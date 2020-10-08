from cognation.formulas.parent import ParentFormula
from cognation.formulas.grandparent import GrandParentFormula
from cognation.formulas.sibling import SiblingFormula
from cognation.formulas.uncle import UncleFormula
from cognation.formulas.cousin import CousinFormula
from cognation.formulas.base import UnknownFormulaException
from cognation.formulas.stepbrother import StepbrotherFormula

FORMULA_TYPE_PARENT = 1
FORMULA_TYPE_GRANDPARENT = 2
FORMULA_TYPE_SIBLING = 3

FORMULA_TYPE_UNCLE = 4
FORMULA_TYPE_COUSIN = 5
FORMULA_TYPE_STEPBROTHER = 6


#  Tells python what to do in every formula case (FORMULA_...)
def formula_builder(type, data):
    normalized_type = int(type)
    if normalized_type == FORMULA_TYPE_PARENT:
        return ParentFormula(data)
    elif normalized_type == FORMULA_TYPE_GRANDPARENT:
        return GrandParentFormula(data)
    elif normalized_type == FORMULA_TYPE_SIBLING:
        return SiblingFormula(data)
    elif normalized_type == FORMULA_TYPE_UNCLE:
        return UncleFormula(data)
    elif normalized_type == FORMULA_TYPE_COUSIN:
        return CousinFormula(data)
    elif normalized_type == FORMULA_TYPE_STEPBROTHER:
        return StepbrotherFormula(data)
    else:
        raise UnknownFormulaException(type)
