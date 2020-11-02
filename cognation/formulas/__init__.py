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

# If in calculation used Hardy-Weinberg law, return 0; else (if used IBD indices) return 1
numToFormula = {
    1: [ParentFormula, 1],
    2: [GrandParentFormula, 1],
    3: [UncleFormula, 1],
    4: [CousinFormula, 1],
    5: [StepbrotherFormula, 1],
    6: [BrotherFormula, 1],
    7: [SiblingFormula, 1],
    8: [TwoParentsFormula, 0],
    9: [TwoChildrenFormula, 0],
    10: [ThreeChildrenFormula, 0]
}


def formula_builder(data_type, data):
    normalized_type = int(data_type)

    if normalized_type in numToFormula.keys():
        target_list = numToFormula[normalized_type]
        return target_list[0](data), target_list[1]
    else:
        raise UnknownFormulaException(data_type)
