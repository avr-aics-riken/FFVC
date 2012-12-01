#ifndef MYGROUP_FACTORY_H
#define MYGROUP_FACTORY_H
#include "common/BBox.h"
#include "polygons/VTree.h"
#include "polygons/Polygons.h"

#include "groups/PolygonGroup.h"
#include "groups/PolygonGroupFactory.h"
#include "CarGroup.h"
#include "Polylib.h"
//#include "BladeGroup.h"
#include <string>
using namespace std;
using namespace PolylibNS;

//class PolygonGroup;

class MyGroupFactory:public PolygonGroupFactory{
 public:
 PolygonGroup* create_instance(std::string class_name){
    if(class_name == CarGroup::get_class_name()){
      //  cout << "CarGroup() called!" <<endl;
      return new CarGroup();
    }

    //  if(class_name == BladeGroup::get_class_name()){
    //    return new BladeGroup();
    //  }

    //default
    return PolygonGroupFactory::create_instance(class_name);

  }
};
#endif // MYGROUPFACTORY_H
