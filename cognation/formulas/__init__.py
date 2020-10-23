from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.cousin import CousinFormula
from cognation.formulas.stepbrother import StepbrotherFormula
from cognation.formulas.brother import BrotherFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.two_parents import TwoParentsFormula
from cognation.formulas.base import UnknownFormulaException


numToFormula = {
    1: ParentFormula,
    2: GrandParentFormula,
    3: UncleFormula,
    4: CousinFormula,
    5: StepbrotherFormula,
    6: BrotherFormula,
    7: SiblingFormula,
    8: TwoParentsFormula
}


def formula_builder(data_type, data):
    normalized_type = int(data_type)

    if normalized_type in numToFormula.keys():
        return numToFormula[normalized_type](data)
    else:
        raise UnknownFormulaException(data_type)
