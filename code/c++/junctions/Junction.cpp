#include "Junction.hpp"

std::ostream& operator<<(std::ostream &os, const Junction &junc) {
  os << junc.p_from->name << "--[" << junc.type() << "]-->" << junc.p_to->name;
  return os;
}
