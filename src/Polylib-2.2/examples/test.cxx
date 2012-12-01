#include "Polylib.h"

using namespace std;
using namespace PolylibNS;


int main(){

  Polylib* pl_instance = Polylib::get_instance();

  pl_instance->load();
  pl_instance->show_group_hierarchy();
  //  pl_instance->show_group_info("car"); // not working??
  string fname="";
  string stl="stl_a";
  string extend="";
  pl_instance->save(&fname,stl,extend);


  
  return 0;


}
