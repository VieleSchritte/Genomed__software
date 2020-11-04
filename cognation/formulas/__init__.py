from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.cousin import CousinFormula
from cognation.formulas.stepbrother import StepbrotherFormula
from cognation.formulas.brother import BrotherFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.two_parents import TwoParentsFormula
from cognation.formulas.base import UnknownFormulaException
from cognation.formulas.two_children import TwoChildrenFormula
from cognation.formulas.three_children import ThreeChildrenFormula
from cognation.formulas.grand_parent_child_yes import ParentGrandChildYes

# If in calculation used Hardy-Weinberg law, return 0; else (if used IBD indices) return 1
numToFormula = {
    1: [ParentFormula, 1],
    2: [TwoChildrenFormula, 0],
    3: [SiblingFormula, 1],
    4: [ThreeChildrenFormula, 0],
    7: [TwoParentsFormula, 0],
    11: [StepbrotherFormula, 1],
    12: [BrotherFormula, 1],
    14: [CousinFormula, 1],
    15: [UncleFormula, 1],
    17: [GrandParentFormula, 1],
    19: [ParentGrandChildYes, 0]
}


def formula_builder(data_type, data):
    normalized_type = int(data_type)

    if normalized_type in numToFormula.keys():
        target_list = numToFormula[normalized_type]
        return target_list[0](data), target_list[1]
    else:
        raise UnknownFormulaException(data_type)
