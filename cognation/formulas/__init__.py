from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.cousin import CousinFormula
from cognation.formulas.stepbrother import StepbrotherFormula
from cognation.formulas.brother import BrotherFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.child_couple import CoupleFormula
from cognation.formulas.base import UnknownFormulaException
from cognation.formulas.two_children import TwoChildrenFormula
from cognation.formulas.three_children import ThreeChildrenFormula
from cognation.formulas.grand_parent_child_yes import ParentGrandChildYes
from cognation.formulas.one_known_supposed import OneKnownSupposedFormula
from cognation.formulas.two_known_supposed import TwoKnownSupposedFormula
from cognation.formulas.two_brothers import TwoBrothersFormula
from cognation.formulas.three_known_supposed import ThreeKnownSupposed
from cognation.formulas.two_couple import TwoCoupleFormula


# If in calculation used Hardy-Weinberg law, return 0; else (if used IBD indices) return 1
numToFormula = {
    1: [ParentFormula, 1],
    2: [TwoChildrenFormula, 0],
    3: [SiblingFormula, 1],
    4: [ThreeChildrenFormula, 0],
    5: [TwoKnownSupposedFormula, 0],
    6: [ThreeKnownSupposed, 0],
    7: [CoupleFormula, 0],
    8: [OneKnownSupposedFormula, 0],
    9: [TwoCoupleFormula, 0],
    10: [],
    11: [StepbrotherFormula, 1],
    12: [BrotherFormula, 1],
    13: [TwoBrothersFormula, 1],
    14: [CousinFormula, 1],
    15: [UncleFormula, 1],
    16: [],
    17: [GrandParentFormula, 1],
    18: [],
    19: [ParentGrandChildYes, 0],
    20: [],
    21: []
}


def formula_builder(data_type, data):
    normalized_type = int(data_type)

    if normalized_type in numToFormula.keys():
        target_list = numToFormula[normalized_type]
        return target_list[0](data), target_list[1]
    else:
        raise UnknownFormulaException(data_type)
