"""Entanglement conformation: routing chains so they actually interlock.

:mod:`~topon.conformation.entanglement.braid` builds a *waypoint braid* --
both partners follow anti-phase ellipses about one shared axis, which makes
the winding count a prescribed input rather than an emergent property of a
tuned bulge.
"""

from topon.conformation.entanglement.allocation import (
    AllocatedContact,
    Allocation,
    ContactRequest,
    Rejection,
    allocate_contacts,
    compose_chain_path,
)
from topon.conformation.entanglement.braid import (
    BraidShape,
    Contact,
    braid_pair,
    braid_path,
    chord_closed_linking,
    far_closed_linking,
    closest_approach,
    feasible_window,
    gap_at,
    linking_number,
    make_contact,
    min_separation,
    plan_braid,
)

__all__ = [
    "AllocatedContact",
    "Allocation",
    "BraidShape",
    "Contact",
    "ContactRequest",
    "Rejection",
    "allocate_contacts",
    "compose_chain_path",
    "braid_pair",
    "braid_path",
    "chord_closed_linking",
    "far_closed_linking",
    "closest_approach",
    "feasible_window",
    "gap_at",
    "linking_number",
    "make_contact",
    "min_separation",
    "plan_braid",
]
