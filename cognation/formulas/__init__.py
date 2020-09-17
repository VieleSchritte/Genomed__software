from .parent import ParentFormula
from .grandparent import GrandParentFormula
from .sibling import SiblingFormula
from .uncle import UncleFormula
from .cousin import CousinFormula
from .stepbrother import StepbrotherFormula
from .base import UnknownFormulaException

FORMULA_TYPE_PARENT = 1
FORMULA_TYPE_GRANDPARENT = 2
FORMULA_TYPE_SIBLING = 3

FORMULA_TYPE_UNCLE = 4
FORMULA_TYPE_COUSIN = 5
FORMULA_TYPE_STEPBROTHER = 6

#Tells python what to do in every formula case (FORMULA_...)
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