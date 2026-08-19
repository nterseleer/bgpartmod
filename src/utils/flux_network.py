"""Currency-parametrised annual flux network of a simulation, as {(src, tgt): flux}.

The generic engine behind the REF-vs-NO-TEP Sankey/flux tables. It builds the same
directed flux network `simulation_manager.carbon_flux_links` produced by hand for carbon,
but for ANY currency (C, N, P, Si) from a single declarative routing table -- so the
per-currency topology irregularities of the model (see below) live in ONE place instead of
being re-hard-coded four times.

Why a table and not pure config introspection: the receiver-side couplings in the config
(`coupled_<proc>_sources`, `coupled_consumers`, ...) give WHO feeds WHOM, but the *currency
routing* lives inside the component code (`if self.name == ...` branches), and it is NOT a
relabelling of the carbon graph. Concretely:
  - lysis      : C -> DOCL, N/P -> DOCS (DON/DOP), Si -> DSi           (three targets!)
  - exudation  : C split DOCS/DOCL, N/P -> DOCS, Si -> DSi
  - DOC->TEP   : C only (TEP carries no N/P/Si)
  - remineralisation / uptake : nutrients only (no detrital C remin; DIC is external)
  - nutrient-only processes: nitrification (NH4->NO3), riverine loads
`_RULES` below is that table; `flux_links` walks it.

Values are read from the saved per-process diagnostic columns
(`<pool>_<sink|source>_<proc>.<cur>`, integrated over the year), so each edge is exact
wherever its term is saved. The two remaining approximations are documented at their
handlers: the grazing split among consumers (`_add_grazing`, preference-weighted until the
per-prey ingestion mirror is saved) and the heterotroph unassimilated split (`_add_hetero
_outputs`, closed by mass balance until `source_ing_*_unassimilated_to_{dom,dim}` are saved).
"""
from typing import Dict, Tuple, List

import numpy as np

from src.utils.simulation_manager import integrate_annual_flux


# Inorganic ("return to DIM") target per currency: respiration -> DIC (C),
# remineralisation -> NH4/DIP/DSi (N/P/Si).
_INORG = {'C': 'DIC', 'N': 'NH4', 'P': 'DIP', 'Si': 'DSi'}

# --- The routing table: one row per direct-edge process family --------------------------
# Each row: (proc-label, emit-role, saved emit-side term, target, currencies).
#   emit-role : 'phy' | 'het' | 'doc' | 'tep' | 'det' | combos ('phy+het', 'dom+det',
#               'all_org') resolved by _emit_pools(); or an explicit list of pool names.
#   term      : the emitting pool's saved term; edge value = integral of <emit>_<term>.<cur>.
#   target    : a fixed node name, or a {cur: node} dict, or _INORG (per-currency DIM node).
#   cur       : tuple of currencies this family carries.
# Exudation, grazing and the heterotroph unassimilated split are NOT here -- they are
# irregular (receiver-side read / consumer split / mass-balance closure) and handled by
# dedicated helpers below.
_ALL = ('C', 'N', 'P', 'Si')
_RULES: List[dict] = [
    dict(proc='respiration', emit='phy',     term='sink_respiration',      tgt=_INORG, cur=('C',)),
    dict(proc='lysis',       emit='phy+het', term='sink_lysis',
         tgt={'C': 'DOCL', 'N': 'DOCS', 'P': 'DOCS', 'Si': 'DSi'},         cur=_ALL),
    dict(proc='mortality',   emit='phy+het', term='sink_mortality',        tgt='DetS', cur=_ALL),
    dict(proc='agg_phy',     emit='phy',     term='sink_aggregation',      tgt='DetL', cur=_ALL),
    dict(proc='agg_doc',     emit='doc',     term='sink_aggregation',      tgt='TEPC', cur=('C',)),
    dict(proc='agg_tep',     emit='tep',     term='sink_aggregation',      tgt='DetL', cur=('C',)),
    dict(proc='agg_dets',    emit=['DetS'],  term='sink_aggregation',      tgt='DetL', cur=_ALL),
    dict(proc='breakdown',   emit='tep',     term='sink_breakdown',        tgt='DOCL', cur=('C',)),
    dict(proc='remin',       emit='dom+det', term='sink_remineralization', tgt=_INORG, cur=('N', 'P', 'Si')),
    dict(proc='vloss',       emit='all_org', term='sink_vertical_loss',    tgt='Export', cur=_ALL),
    dict(proc='leak',        emit='all_org', term='sink_leakage_out',      tgt='Leak', cur=_ALL),
]


def _roles(config) -> Dict[str, List[str]]:
    """Group the organic/inorganic state variables by model role, read off the config
    classes (so pool NAMES are not hard-coded, only the special TEP split by name)."""
    r = {'phy': [], 'het': [], 'doc': [], 'tep': [], 'det': [], 'dim': []}
    for key, cfg in config.items():
        if key == 'formulation':
            continue
        cn = cfg['class'].__name__
        if cn == 'Phyto':
            r['phy'].append(key)
        elif cn == 'Heterotrophs':
            r['het'].append(key)
        elif cn == 'Detritus':
            r['det'].append(key)
        elif cn == 'DIM':
            r['dim'].append(key)
        elif cn == 'DOM':
            (r['tep'] if key == 'TEPC' else r['doc']).append(key)
    return r


def _emit_pools(spec, roles) -> List[str]:
    """Resolve a rule's emit-role spec into the list of emitting pool names."""
    if isinstance(spec, (list, tuple)):
        return list(spec)
    combos = {
        'phy+het': roles['phy'] + roles['het'],
        'dom+det': roles['doc'] + roles['tep'] + roles['det'],
        'all_org': roles['phy'] + roles['het'] + roles['doc'] + roles['tep'] + roles['det'],
    }
    return combos.get(spec, roles.get(spec, []))


def _tgt(tgt, cur):
    """Resolve a rule target for one currency (fixed node or per-currency dict)."""
    return tgt[cur] if isinstance(tgt, dict) else tgt


def _graze_split(prey, roles, config):
    """Static-preference split of a prey's grazing loss among its consumers (heterotrophs).

    The realised preference is biomass-weighted (Heterotrophs.get_source_ingestion), so this
    is the ONE approximation in the network -- exact once the per-prey ingestion mirror
    (`<consumer>_ing_<prey>_C`) is saved; see module docstring."""
    w = {h: config[h]['parameters'].get(f'pref_{prey}', 0.0) for h in roles['het']}
    w = {h: p for h, p in w.items() if p > 0}
    tot = sum(w.values())
    return {h: p / tot for h, p in w.items()} if tot > 0 else {}


def flux_links(sim, currency: str = 'C', period: int = 2023) -> Dict[Tuple[str, str], float]:
    """Depth-integrated annual flux network {(src, tgt): flux} for one currency.

    Args:
        sim: a Model simulation (needs .df with the saved process columns and .config).
        currency: 'C' (default), 'N', 'P' or 'Si'. The C path reproduces
                  `simulation_manager.carbon_flux_links` exactly (regression-tested);
                  N/P/Si follow the same table but require the element-resolved diagnostics
                  (a dedicated re-run) to be fully populated.
        period: year to integrate over.

    Returns:
        {(src, tgt): annual flux [mmol <cur> m-2 yr-1]}; external nodes: DIC/NH4/DIP/DSi
        (return to inorganic), Export (settling), Leak (calibrated kleak loss).
    """
    roles = _roles(sim.config)
    cols = set(sim.df.columns)
    A = lambda col: (integrate_annual_flux(sim, col, period) if col in cols else 0.0)
    L: Dict[Tuple[str, str], float] = {}

    def add(src, tgt, val):
        if val and not (isinstance(val, float) and np.isnan(val)):
            L[(src, tgt)] = L.get((src, tgt), 0.0) + val

    # 1) Direct-edge rules (value read from the emitting pool's saved term).
    for rule in _RULES:
        if currency not in rule['cur']:
            continue
        for src in _emit_pools(rule['emit'], roles):
            add(src, _tgt(rule['tgt'], currency), A(f"{src}_{rule['term']}.{currency}"))

    # 2) Exudation (irregular: read receiver-side, C splits DOCS/DOCL).
    _add_exudation(add, A, roles, currency)

    # 3) Grazing (prey-side exact total, split among consumers exactly from the ingestion
    #    mirror when available, else by static preference).
    _add_grazing(add, A, cols, roles, sim.config, currency)

    # 4) Heterotroph unassimilated outputs, closed by mass balance from the routed intake.
    _add_hetero_outputs(add, A, L, roles, sim.config, currency)

    return L


def _add_exudation(add, A, roles, currency):
    """Phytoplankton exudation. C: split small/large read from DOCS/DOCL source side
    (the split fraction is internal to Phy); N/P: all to DOCS; Si: to DSi (emit-side, so it
    does not double-count the lysis-Si already routed to DSi by the lysis rule)."""
    for phy in roles['phy']:
        if currency == 'C':
            add(phy, 'DOCS', A('DOCS_source_exudation.C'))
            add(phy, 'DOCL', A('DOCL_source_exudation.C'))
        elif currency in ('N', 'P'):
            add(phy, 'DOCS', A(f'DOCS_source_exudation.{currency}'))
        elif currency == 'Si':
            add(phy, 'DSi', A(f'{phy}_sink_exudation.Si'))


def _add_grazing(add, A, cols, roles, config, currency):
    """Grazing: exact on the prey side (`<prey>_sink_ingestion.<cur>`); the split among
    consumers is EXACT when the per-prey ingestion mirror `<cons>_ing_<prey>_C` is saved
    (each consumer's realised C-share, which the N/P/Si totals also follow), and falls back
    to the static-preference split otherwise (older sims without the mirror)."""
    prey_pools = roles['doc'] + roles['tep'] + roles['det'] + roles['het']
    for prey in prey_pools:
        tot = A(f'{prey}_sink_ingestion.{currency}')
        if not tot:
            continue
        cons = [h for h in roles['het'] if config[h]['parameters'].get(f'pref_{prey}', 0.0) > 0]
        mir_cols = [f'{c}_ing_{prey}_C' for c in cons]
        if cons and all(mc in cols for mc in mir_cols):        # exact realised split
            mir = {c: A(mc) for c, mc in zip(cons, mir_cols)}
            s = sum(mir.values())
            shares = {c: mir[c] / s for c in cons} if s > 0 else {}
        else:                                                  # fallback: static preference
            shares = _graze_split(prey, roles, config)
        for c, w in shares.items():
            add(prey, c, tot * w)


def _add_hetero_outputs(add, A, L, roles, config, currency):
    """Heterotroph metabolic outputs, closed by mass balance from the routed intake.

    C/N/P: the unassimilated fraction splits to DOCS (sloppy feeding, f_unass_excr) and to
    the inorganic pool (sloppy-to-DIM + respiration for C). Si is different -- heterotrophs
    carry no Si pool, so ALL ingested Si is unassimilated and splits to DetS (f_unass_Si) and
    DSi (1 - f_unass_Si), never to DOCS. Lysis/mortality/vertical loss are emitted by the
    direct rules."""
    for h in roles['het']:
        intake = sum(v for (a, b), v in L.items() if b == h)  # grazing into h (this currency)
        p = config[h]['parameters']
        if currency == 'Si':
            f_si = p.get('f_unass_Si', 0.9)
            add(h, 'DetS', intake * f_si)
            add(h, 'DSi', intake * (1.0 - f_si))
        else:
            assim = A(f'{h}_source_ing_{currency}_assimilated')
            resp = A(f'{h}_sink_respiration.{currency}') if currency == 'C' else 0.0
            f_dom = p['f_unass_excr']
            unassim = max(intake - assim - resp, 0.0)
            add(h, 'DOCS', unassim * f_dom)
            add(h, _INORG[currency], unassim * (1 - f_dom) + resp)
