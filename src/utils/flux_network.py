"""Currency-parametrised annual flux network of a simulation, as {(src, tgt): flux}.

The generic engine behind the REF-vs-NO-TEP Sankey/flux tables. It builds the same
directed flux network `carbon_flux_links` (below) produced by hand for carbon,
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
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd

from src.config_model import varinfos
from src.core.legacy import NON_COMPONENT_KEYS


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
        if key in NON_COMPONENT_KEYS:
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
                  `carbon_flux_links` (below) exactly (regression-tested);
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


# ============================================================================
# ANNUAL FLUX INTEGRATION AND COMPARISON
# ============================================================================
# Reading side of the network above: integrate a saved flux term over a year, and
# compare two simulations term by term.

def integrate_annual_flux(simulation, var, period: int = 2023) -> float:
    """Depth-integrated annual sum of a rate flux -> [<flux unit> m-2 yr-1].

    Generalises _integrate_PP to any rate diagnostic (or a sum of several). Fluxes in the
    model are volumetric rates [X m-3 d-1]; multiplying the annual time-integral by the
    water-column depth gives an areal annual budget directly comparable across pathways
    and across simulations (the currency of a carbon-flux Sankey / bar comparison).

    Args:
        simulation: Model with .df (DatetimeIndex) and .setup.base_water_depth.
        var: Single df column name, or a list of column names whose values are summed
             (e.g. the four aggregate-coupled sink_vertical_loss.C for total C export).
        period: Year to integrate over. If None, integrate the whole run.

    Returns:
        Depth-integrated annual flux [<flux unit> m-2 yr-1]. Missing columns contribute
        NaN (so a partially-diagnosed run surfaces the gap rather than silently dropping
        a term).
    """
    df = simulation.df[simulation.df.index.year == period] if period else simulation.df
    dt = (df.index[1] - df.index[0]).total_seconds() / 86400  # timestep in days
    depth = simulation.setup.base_water_depth                  # water column depth in m
    cols = [var] if isinstance(var, str) else list(var)
    total = 0.0
    for c in cols:
        if c not in df.columns:
            return float('nan')
        total += pd.to_numeric(df[c], errors='coerce').sum()
    return total * dt * depth


def compare_annual_fluxes(sims, fluxes, period: int = 2023, names=None,
                          ref_index: int = 0) -> pd.DataFrame:
    """Tabulate depth-integrated annual C fluxes for several simulations and their change.

    The Step-1 signal check behind the REF vs NO-TEP carbon-cycle analysis: does switching
    the TEP->flocculation coupling off reorganise the fluxes (PP -> DOC/Det -> bacteria ->
    export), not just the standing stocks? Feed it the pathway fluxes and read the relative
    change column before deciding whether a Sankey is worth building.

    Args:
        sims: list of Model simulations (e.g. [sim_ref, sim_notep]).
        fluxes: either a list of df column names, or a dict {label: column-or-list} where a
                list value is summed (integrate_annual_flux). Labels keep the table readable.
        period: year to integrate over (default 2023), passed to integrate_annual_flux.
        names: display names per simulation (default: each sim.name).
        ref_index: which simulation is the reference for the relative-change column.

    Returns:
        DataFrame indexed by flux label, one column of [mmol C m-2 yr-1] per simulation
        plus 'rel_change_%' ((other - ref)/ref * 100 for the two-simulation case, else the
        last simulation vs the reference).
    """
    if names is None:
        names = [getattr(s, 'name', f'sim{i}') for i, s in enumerate(sims)]
    if not isinstance(fluxes, dict):
        fluxes = {v if isinstance(v, str) else '+'.join(v): v for v in fluxes}

    data = {name: [integrate_annual_flux(s, var, period) for var in fluxes.values()]
            for name, s in zip(names, sims)}
    table = pd.DataFrame(data, index=list(fluxes.keys()))

    ref_col = names[ref_index]
    other_col = names[-1] if names[-1] != ref_col else names[min(1, len(names) - 1)]
    with np.errstate(divide='ignore', invalid='ignore'):
        table['rel_change_%'] = (table[other_col] / table[ref_col] - 1.0) * 100.0
    return table


def _integrate_PP(simulation, var: str = 'Phy_source_PP.C', period: int = 2023) -> float:
    """Integrate primary production over a given year.

    Args:
        simulation: Model simulation object with .df (DatetimeIndex) and .setup attributes.
        var: Name of the PP rate variable in simulation.df [mmol C m-3 d-1].
        period: Year to integrate over. If None, uses the full simulation.

    Returns:
        Depth-integrated annual PP [mmol C m-2].
    """
    return integrate_annual_flux(simulation, var, period)


def compare_annual_PP(sim1, sim2, period: int = 2023, var: str = 'Phy_source_PP.C') -> None:
    """Print depth-integrated annual PP for two simulations and their relative difference.

    Args:
        sim1: Reference simulation.
        sim2: Simulation to compare against sim1.
        period: Year to integrate over (default: 2023).
        var: PP rate variable name in simulation.df [mmol C m-3 d-1].
    """
    pp1 = _integrate_PP(sim1, var, period)
    pp2 = _integrate_PP(sim2, var, period)
    rel_change = (pp2 / pp1 - 1) * 100
    direction = 'increase' if rel_change > 0 else 'decrease'
    var_label = varinfos.doutput.get(var, {}).get('longname', var)
    print(f"{var_label} ({period}):  {sim1.name} = {pp1:.1f}  |  {sim2.name} = {pp2:.1f}  [mmol C m-2 yr-1]"
          f"  ->  {rel_change:+.1f}% {direction}")


# Canonical left->right, top->bottom reading order of the flux network (matches the Sankey
# layout), incl. the aggregate node names (DOC/Det/Bac) so aggregated tables order too.
_SANKEY_POOL_ORDER = ['Phy', 'DOCS', 'DOCL', 'DOC', 'TEPC', 'DetS', 'DetL', 'Det',
                      'BacF', 'BacA', 'Bac', 'HF', 'Cil']


def carbon_flux_links(sim, period: int = 2023) -> Dict[Tuple[str, str], float]:
    """Depth-integrated annual carbon-flux network of one simulation, as {(src, tgt): flux}.

    Thin currency='C' wrapper around `flux_network.flux_links` (kept for backward
    compatibility and as the regression anchor: the generic engine reproduces this network
    edge-for-edge). Every link is a depth-integrated annual C flux [mmol C m-2 yr-1] between
    two organic-C pools (Phy, DOCS, DOCL, TEPC, DetS, DetL, BacF, BacA, HF, Cil) or an
    external node: 'DIC' (respiration + remineralization + sloppy-feeding-to-DIM), 'Export'
    (settling sink_vertical_loss) and 'Leak' (the calibrated kleak loss). Primary production
    is not drawn as an inflow, so Phy is the network source and its bar height reads as PP.

    Topology, the grazing preference split (the only approximation) and the heterotroph
    mass-balance closure are documented in `flux_network`. Node totals reconcile with the
    diagnosed C_sources/C_sinks; the system closes to ~2-3% of PP (grazing-split residual).
    """
    return flux_links(sim, currency='C', period=period)


def aggregate_flux_links(links, groups) -> Dict[Tuple[str, str], float]:
    """Merge nodes of a flux network {(src, tgt): flux} into aggregate nodes.

    Companion to `flux_balance_table`/`carbon_flux_links` for a coarser view (e.g. DOC =
    DOCS+DOCL, Det = DetS+DetL, Bac = BacF+BacA). Parallel links that collapse onto the same
    (src, tgt) after relabelling are summed; links internal to a single aggregate become
    self-loops and are dropped (they no longer cross a node boundary, so they are not a flux
    between the reported pools -- e.g. DOCS->DOCL vanishes inside 'DOC').

    Args:
        links: the flux network to coarsen.
        groups: {member_node: aggregate_name}; nodes absent from it keep their own name.

    Returns:
        A new links dict on the aggregated node set.
    """
    m = lambda n: groups.get(n, n)
    out: Dict[Tuple[str, str], float] = {}
    for (a, b), v in links.items():
        A, B = m(a), m(b)
        if A == B:
            continue
        out[(A, B)] = out.get((A, B), 0.0) + v
    return out


def flux_balance_table(links_a, links_b, names=('REF', 'NO-TEP'), pools=None,
                       min_flux=0.0) -> pd.DataFrame:
    """Per-pool incoming/outgoing flux balance of two flux-link networks + their change.

    The tabular companion to the REF-vs-NO-TEP Sankey: given two link dicts
    {(src, tgt): flux} as returned by `carbon_flux_links` (currency-agnostic -- works for a
    future N/P/DSi network too), it lists, for each pool, every flux entering and leaving it,
    with the absolute value, its share of that pool's total in- (or out-) throughput, in both
    simulations, and the absolute + relative change between them. Reads straight off the
    Sankey links so the table and the diagram tell exactly the same story.

    Args:
        links_a, links_b: the two networks; `names[0]` is the reference for the change columns.
        names: (name_a, name_b) column labels.
        pools: pool order to report. If None, every non-external node that has any flux, in
               the canonical Sankey reading order (`_SANKEY_POOL_ORDER`; unknown nodes appended
               by descending throughput). External nodes (DIC, Export, Leak, DIC_in, CO2) only
               ever appear as counterparts, never as a pool.
        min_flux: drop counterpart links whose value is < min_flux in BOTH sims (the per-pool
                  '(total)' rows always use the full, unfiltered throughput).

    Returns:
        Tidy DataFrame, columns:
          pool, direction ('in'|'out'), counterpart,
          <name_a>, <name_a>_pct, <name_b>, <name_b>_pct, abs_change, rel_change_%
        A '(total)' counterpart row precedes each (pool, direction) block (pct = 100),
        giving the pool's absolute in/out throughput and how much it moved.
    """
    na, nb = names
    _EXTERNAL = {'DIC', 'CO2', 'Export', 'Leak', 'DIC_in'}
    nodes = set(a for a, _ in links_a) | set(b for _, b in links_a) \
        | set(a for a, _ in links_b) | set(b for _, b in links_b)

    def out_tot(L, p):
        return sum(v for (a, b), v in L.items() if a == p)

    def in_tot(L, p):
        return sum(v for (a, b), v in L.items() if b == p)

    if pools is None:
        cand = [p for p in nodes if p not in _EXTERNAL]
        rank = {p: i for i, p in enumerate(_SANKEY_POOL_ORDER)}
        pools = sorted(cand, key=lambda p: (rank.get(p, len(rank)),
                                            -(out_tot(links_a, p) + in_tot(links_a, p))))

    rows = []
    for pool in pools:
        for direction in ('in', 'out'):
            if direction == 'out':
                tot_a, tot_b = out_tot(links_a, pool), out_tot(links_b, pool)
                counterparts = {b for (a, b) in links_a if a == pool} \
                    | {b for (a, b) in links_b if a == pool}
                val = lambda L, cp: L.get((pool, cp), 0.0)
            else:
                tot_a, tot_b = in_tot(links_a, pool), in_tot(links_b, pool)
                counterparts = {a for (a, b) in links_a if b == pool} \
                    | {a for (a, b) in links_b if b == pool}
                val = lambda L, cp: L.get((cp, pool), 0.0)
            if tot_a == 0.0 and tot_b == 0.0:
                continue

            def mkrow(cp, va, vb, ta, tb):
                return {
                    'pool': pool, 'direction': direction, 'counterpart': cp,
                    na: va, f'{na}_pct': 100.0 * va / ta if ta else np.nan,
                    nb: vb, f'{nb}_pct': 100.0 * vb / tb if tb else np.nan,
                    'abs_change': vb - va,
                    'rel_change_%': (vb / va - 1.0) * 100.0 if va else np.nan,
                }

            rows.append(mkrow('(total)', tot_a, tot_b, tot_a, tot_b))
            cps = sorted(counterparts, key=lambda cp: -val(links_a, cp))
            for cp in cps:
                va, vb = val(links_a, cp), val(links_b, cp)
                if max(va, vb) < min_flux:
                    continue
                rows.append(mkrow(cp, va, vb, tot_a, tot_b))
    return pd.DataFrame(rows)
