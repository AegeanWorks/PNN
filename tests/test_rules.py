"""Unit tests for the property rules (PRINCIPLE-014: failing-first boundaries)."""

import unittest

from pnn.compounds import Compound
from pnn.rules import FAIL, PASS, UNKNOWN, evaluate


def ideal(**overrides) -> Compound:
    values = {
        "compound_id": "CMP-001",
        "mw": 350.0,
        "logp": 2.5,
        "tpsa": 55.0,
        "hbd": 2.0,
        "hba": 5.0,
        "rotatable_bonds": 4.0,
        "charge": "neutral",
        "pka": 8.0,
        "logbb": -0.2,
        "brain_plasma_ratio": 0.5,
        "oral_bioavailability": 0.35,
        "cyp3a4_inhibition": 10.0,
        "herg_liability": 5.0,
        "pnn_ec50_nm": 5.0,
    }
    values.update(overrides)
    return Compound(**values)


def by_name(compound: Compound) -> dict:
    return {criterion.name: criterion for criterion in evaluate(compound)}


class PropertyRuleTest(unittest.TestCase):
    def test_ideal_compound_passes_every_rule(self) -> None:
        criteria = by_name(ideal())
        self.assertTrue(all(item.status == PASS for item in criteria.values()))
        self.assertTrue(all(item.score == 1.0 for item in criteria.values()))

    def test_molecular_weight_boundaries(self) -> None:
        self.assertEqual(by_name(ideal(mw=200.0))["mw"].score, 0.7)
        self.assertEqual(by_name(ideal(mw=300.0))["mw"].score, 1.0)
        self.assertEqual(by_name(ideal(mw=450.0))["mw"].score, 0.7)
        self.assertEqual(by_name(ideal(mw=451.0))["mw"].status, FAIL)

    def test_tpsa_and_rotor_ceilings(self) -> None:
        criteria = by_name(ideal(tpsa=65.0, rotatable_bonds=8.0))
        self.assertEqual(criteria["tpsa"].score, 0.7)
        self.assertEqual(criteria["rotatable_bonds"].status, FAIL)

    def test_missing_optional_properties_are_unknown(self) -> None:
        criteria = by_name(
            Compound(
                compound_id="CMP-002",
                mw=350.0,
                logp=2.5,
                tpsa=55.0,
                hbd=2.0,
                hba=5.0,
                rotatable_bonds=4.0,
            )
        )
        self.assertEqual(criteria["logbb"].status, UNKNOWN)
        self.assertEqual(criteria["herg_liability"].status, UNKNOWN)
        self.assertEqual(criteria["pka"].status, UNKNOWN)

    def test_basic_pka_outside_window_fails(self) -> None:
        criteria = by_name(ideal(pka=6.0, charge="basic"))
        self.assertEqual(criteria["pka"].status, FAIL)

    def test_unrecognized_charge_is_unknown(self) -> None:
        criteria = by_name(ideal(charge="cationic-ish"))
        self.assertEqual(criteria["charge"].status, UNKNOWN)

    def test_oral_bioavailability_accepts_fraction_and_percent(self) -> None:
        self.assertAlmostEqual(by_name(ideal(oral_bioavailability=0.25))["oral_bioavailability"].value, 25.0)
        self.assertAlmostEqual(by_name(ideal(oral_bioavailability=25.0))["oral_bioavailability"].value, 25.0)
        self.assertEqual(by_name(ideal(oral_bioavailability=0.05))["oral_bioavailability"].status, FAIL)

    def test_hard_liabilities_fail(self) -> None:
        criteria = by_name(ideal(herg_liability=45.0, cyp3a4_inhibition=70.0, logbb=-1.5))
        self.assertEqual(criteria["herg_liability"].status, FAIL)
        self.assertEqual(criteria["cyp3a4_inhibition"].status, FAIL)
        self.assertEqual(criteria["logbb"].status, FAIL)


if __name__ == "__main__":
    unittest.main()
