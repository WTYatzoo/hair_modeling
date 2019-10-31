#include "simulator.h"
#include "io.h"
using namespace std;

simulator::simulator(boost::property_tree::ptree &para_tree)
{
  kind=para_tree.get<string>("kind.value");
  if(kind=="default")
    {
      default_simulate(para_tree);
    }
  else if(kind=="use_data_from_file")
    {
      simulate(para_tree);
    }
}

int simulator::simulate(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  out_dir=para_tree.get<string>("out_dir.value"); //此处要给到最深一层 // such as ./example/mybeam3/ {mybeam3_fine_constitutive_model_ori mybeam3_coarsen_constitutive_model_ori mybeam3_coarsen_constitutive_model_new}
  
  simulation_type=para_tree.get<string>("simulation.value","static"); // 
  frame=para_tree.get<int>("frame.value",3000); //
  dt=para_tree.get<double >("dt.value");  //
  gravity=para_tree.get<double>("gravity.value"); //
  density=para_tree.get<double>("density.value"); //
  line_search=para_tree.get<int >("line_search.value"); //
  weight_line_search=para_tree.get<double >("weight_line_search.value"); //
  constitutive_model=para_tree.get<string>("constitutive_model.value"); //constitutive_model & input_mat_kind couple

  string object_name=para_tree.get<string>("object_name.value"); //
  string fine_coarsen=para_tree.get<string>("fine_coarsen.value"); // fine or coarsen
  string input_object=para_tree.get<string>("input_object.value"); // 具体到.vtk文件
  string input_mat=para_tree.get<string>("input_mat.value"); // 具体到.mat文件
  string input_constraint=para_tree.get<string>("input_constraint.value"); //具体到.csv文件
  int input_mat_kind=para_tree.get<int>("input_mat_kind.value");// 1: material_para 0:stiffness_tensor // int

  int level=para_tree.get<int>("level.value");
  dmetric=myio.getDmetric(input_object);
  if(fine_coarsen=="fine")
    {
      ;
    }
  else if(fine_coarsen=="coarsen")
    {
      ;
    }
  force_function=para_tree.get<string>("force_function.value");
  
  if(simulation_type=="static")
    {
      //make sure that it is a static simulator 
      dt=100;
      frame=1;
    }
  else if(simulation_type=="dynamic")
    {
      // do not modify the para from bash
    }
  size_t i,j,k;
  double energy;

  object* test_obj=new object(dmetric,input_object,dt,density,line_search,weight_line_search,constitutive_model); 
  myio.getMatPara(test_obj,input_mat,input_mat_kind);
  myio.getConstraintFromCsv(test_obj,input_constraint);

  {
    if(force_function=="gravity")
      {
	
	for(i=0;i<test_obj->num_vertex;++i)
	  {
	    test_obj->myvertexs[i].force_external=myvector(0,-1*gravity*0.98,0)*test_obj->myvertexs[i].mass;
	    //    test_obj->myvertexs[i].force_external=myvector(0,0,-1*gravity)*test_obj->myvertexs[i].mass;
	  }
	  /*
	  test_obj->myvertexs[4].force_external=myvector(1000,0,0)*test_obj->myvertexs[4].mass;
	  test_obj->myvertexs[5].force_external=myvector(1000,0,0)*test_obj->myvertexs[5].mass;
	  test_obj->myvertexs[6].force_external=myvector(1000,0,0)*test_obj->myvertexs[6].mass;
	  test_obj->myvertexs[7].force_external=myvector(1000,0,0)*test_obj->myvertexs[7].mass;*/
	  
	//生成中间往上的重力，两边往下的重力，for mybeam7, try to find one weakness of numerical coarsening 但是因为rectangular function 不连续，导致力并不等价，所以结果不具有意义
	/*
	{
	  double x_center=0;
	  for(i=0;i<test_obj->num_vertex;i++)
	    {
	      x_center+=test_obj->myvertexs[i].location(0);
	    }
	  x_center/=(double)test_obj->num_vertex;
	  for(i=0;i<test_obj->num_vertex;i++)
	    {
	      if(fabs(test_obj->myvertexs[i].location(0)-x_center)<0.75) //因为这个case的单位长度精确到小数点后一位,所以我用个0.75避免因为浮点误差带来的不精确囊括
		{
		  test_obj->myvertexs[i].force_external=myvector(0,0,3*gravity)*test_obj->myvertexs[i].mass;
		}
	      else
		{
		  test_obj->myvertexs[i].force_external=myvector(0,0,-3*gravity)*test_obj->myvertexs[i].mass;
		}
	    }
	    }*/
	//用rational function 代替rectangular function
	/*
	{
	  double x_center=0;
	  double x_min=test_obj->myvertexs[0].location(0);
	  for(i=0;i<test_obj->num_vertex;i++)
	    {
	      if(test_obj->myvertexs[i].location(0)<x_min)
		{
		  x_min=test_obj->myvertexs[i].location(0);
		}
	      x_center+=test_obj->myvertexs[i].location(0);
	    }
	  x_center/=(double)test_obj->num_vertex;
	  double len=fabs(x_center-x_min);
	  for(i=0;i<test_obj->num_vertex;i++)
	    {
	      test_obj->myvertexs[i].force_external=myvector(0,0,gravity)*test_obj->myvertexs[i].mass*(1.0/(pow(4*(test_obj->myvertexs[i].location(0)-x_center)/len,8)+1)-0.5)*6;
	    }
	    }*/
      }
    else if(force_function=="rotation")
      {
	for(i=0;i<test_obj->num_vertex;++i)
	  {
	    myvector help=test_obj->myvertexs[i].location-test_obj->center_loc;
	    myvector torque=myvector(gravity,0,0);
	    myvector r=help-(help.dot(torque))*torque;
	    test_obj->myvertexs[i].force_external=torque.cross(r)*test_obj->myvertexs[i].mass;
	  }
      }
    else if(force_function=="compression")
      {
	/*
	for(i=0;i<test_obj->num_vertex;++i)
	  {
	    myvector help=test_obj->center_loc-test_obj->myvertexs[i].location;
	    double scale=5.0*test_obj->myvertexs[i].mass;
	    test_obj->myvertexs[i].force_external=gravity*scale*help.dot(myvector(1.0,0,0))*myvector(1.0,0,0);
	    cout<<test_obj->myvertexs[i].mass*10000<<endl;
	    }*/
      }
  }

  test_obj->checkFixedOrFree();
  for(i=0;i<frame;++i)
    {
      test_obj->dynamicSimulator();
      energy=test_obj->calElasticEnergy();
      printf("energy is :%lf \n ",energy);
      stringstream ss;string frame_str;ss<<i; ss>>frame_str;
      stringstream ee; string gravity_str;ee<<gravity; ee>>gravity_str;


      if(simulation_type=="dynamic")
	{
	  string path_name=out_dir+"/"+object_name+"_"+force_function+"_"+simulation_type+"_"+gravity_str+"_"+frame_str+".vtk";
	  myio.saveAsVTKwithPara(test_obj,path_name);
	  path_name=out_dir+"/"+object_name+"_"+force_function+"_"+simulation_type+"_"+gravity_str+"_"+frame_str+".csv";
	  myio.saveTimeNormPair(test_obj,path_name);
	}
      else if(simulation_type=="static")
	{
	  string path_name=out_dir+"/"+object_name+"_"+force_function+"_"+simulation_type+"_"+frame_str+"_"+gravity_str+".vtk";
	  myio.saveAsVTKwithPara(test_obj,path_name);
	  path_name=out_dir+"/"+object_name+"_"+force_function+"_"+simulation_type+"_"+frame_str+"_"+gravity_str+".csv";
	  myio.saveTimeNormPair(test_obj,path_name);
	}
    }
  delete test_obj;
  
  return 0;
}

int simulator::default_simulate(boost::property_tree::ptree &para_tree)
{
  io myio=io();
  out_dir=para_tree.get<string>("out_dir.value"); //
  simulation_type=para_tree.get<string>("simulation.value","static"); // 
  frame=para_tree.get<int>("frame.value",3000); //
  dt=para_tree.get<double >("dt.value");  //
  gravity=para_tree.get<double>("gravity.value"); //
  density=para_tree.get<double>("density.value"); //
  line_search=para_tree.get<int >("line_search.value"); //
  weight_line_search=para_tree.get<double >("weight_line_search.value"); //
  constitutive_model=para_tree.get<string>("constitutive_model.value"); //
  force_function="gravity";  //
  dmetric=0.01;//

  double PoissonRatio=para_tree.get<double >("poissonratio.value",0.45); 
  double YoungModulus=para_tree.get<double >("youngmodulus.value",1e5);
  
  if(simulation_type=="static")
    {
      //make sure that it is a static simulator 
      dt=100;
      frame=1;
    }
  else if(simulation_type=="dynamic")
    {
      // do not modify the para from bash
    }
  size_t i,j,k;
  double material[2];
  double energy;
  
  material[0]=YoungModulus/(1+PoissonRatio)*0.5;
  material[1]=YoungModulus*PoissonRatio/(1+PoissonRatio)/(1-2*PoissonRatio);

  object* test_obj=new object(dmetric,12,2,2,dt,density,line_search,weight_line_search,constitutive_model);
  configMaterial(test_obj,material);
  configForce(test_obj);
  test_obj->checkFixedOrFree();
  for(i=0;i<frame;++i)
    {
      test_obj->dynamicSimulator();
      energy=test_obj->calElasticEnergy();
      printf("energy is :%lf \n ",energy);
      stringstream ss;string frame_str;ss<<i; ss>>frame_str;
      string path_name=out_dir+"/test_obj_"+simulation_type+"_"+frame_str+".vtk";
      myio.saveAsVTK(test_obj,path_name);
    }
  delete test_obj;
  
  return 0;
}

int simulator::configMaterial(object *myobject,const double (&material_para)[2])
{
  size_t num_hexahedrons=myobject->num_hexahedrons;
  size_t i,j,k,a,b,c;
  for(i=0;i<num_hexahedrons;++i)
    {
      for(a=0;a<2;++a)
	{
	  for(b=0;b<2;++b)
	    {
	      for(c=0;c<2;++c)
		{
		  myobject->myhexahedrons[i].material_para[a][b][c][0]=material_para[0];
		  myobject->myhexahedrons[i].material_para[a][b][c][1]=material_para[1];
		}
	    }
	}
      myobject->myhexahedrons[i].stiffness_tensor(0,0)=myobject->myhexahedrons[i].stiffness_tensor(1,1)=myobject->myhexahedrons[i].stiffness_tensor(2,2)=material_para[0]*2+material_para[1];
      myobject->myhexahedrons[i].stiffness_tensor(0,1)=myobject->myhexahedrons[i].stiffness_tensor(0,2)=myobject->myhexahedrons[i].stiffness_tensor(1,2)=myobject->myhexahedrons[i].stiffness_tensor(1,0)=myobject->myhexahedrons[i].stiffness_tensor(2,0)=myobject->myhexahedrons[i].stiffness_tensor(2,1)=material_para[1];
      myobject->myhexahedrons[i].stiffness_tensor(3,3)=myobject->myhexahedrons[i].stiffness_tensor(4,4)=myobject->myhexahedrons[i].stiffness_tensor(5,5)=material_para[0];
    }
  return 0;
}

int simulator::configForce(object *myobject)
{
  size_t i,j,k;
  size_t lc=myobject->lc; size_t wc=myobject->wc; size_t hc=myobject->hc;
  
  myobject->num_fixed=0;
  for(j=0;j<wc;++j)
    {
      for(k=0;k<hc;++k)
	{
	  myobject->myvertexs[myobject->index_for_vertex[0][j][k]].isFixed=1;
	  ++myobject->num_fixed;
	}
    }
  if(force_function=="gravity")
    {
      for(i=0;i<lc;++i)
	{
	  for(j=0;j<wc;++j)
	    {
	      for(k=0;k<hc;++k)
		{
		  myobject->myvertexs[myobject->index_for_vertex[i][j][k]].force_external=myvector(0,0,-1*gravity)*myobject->myvertexs[myobject->index_for_vertex[i][j][k]].mass;
		}
	    } 
	}
    }
  return 0;
}
