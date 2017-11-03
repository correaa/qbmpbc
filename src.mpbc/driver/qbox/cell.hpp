#ifndef QBOX_CELL_HPP
#define QBOX_CELL_HPP
#include "UnitCell.h"
namespace qbox{
  //	typedef UnitCell cell;
  struct cell : UnitCell{
    cell(UnitCell const& other) : UnitCell(other){}
    cell(D3vector const& v1, D3vector const& v2, D3vector const& v3) : 
      UnitCell(v1, v2, v3){}
    qbox::volume volume() const{
      return UnitCell::volume()*qbox::volume::unit_type();
    }
    cell& operator*=(double scale){
      D3vector a0=a(0), a1=a(1), a2=a(2);
      this->set(
		a0*scale,
		a1*scale,
		a2*scale
      );
      return *this;
    }
  };
}
#endif
