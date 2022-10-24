#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>


#include "hdf5_save.h"
#define H5_USE_18_API 
#define H5Acreate_vers 2
#define H5Dcreate_vers 2
#define H5Dopen_vers 2


#include <hdf5.h>


void hdf5_save_var(int *sdiv, int *mdiv,
		   char fname[],
		   /* ------------------- */
                   char eos_type[],
                   char eos_file[],
                   double *eos_k,
                   double *Gamma_P,
		   /* ------------------- */
                   char   rotation_type[],
                   double *A_diff,
		   /* ------------------- */
                   double *r_ratio,
                   double *rho0_center,
                   double *r_e,	
		   /* ------------------- */
                   double *s_qp,  /*SDIV+1*/
                   double *mu,	  /*MDIV+1*/
		   /* ------------------- */
                   double **rho_potential, 
                   double **gama, 
                   double **alpha,
                   double **omega,
                   double **energy,
                   double **pressure,
                   double **enthalpy,
                   double **velocity_sq,
                   double *Omega,
                   double *Omega_e,
		   double **Omega_diff)
{

  hid_t       file_id;   /* file identifier */
  herr_t      status;

  hid_t       attribute_id; /* identifiers for attributes*/
  hid_t       attributeH5type,attributeH5c128;

  hid_t       dataset_id, dataspace_id;  /* identifiers for dsets*/
  hsize_t     dims[3];

  int         i,j;
  double      *dset_data;
  double      **var;
  int         varIndex;

  /* =========================================== */
  /* Create a new file using default properties. */
  /* =========================================== */
  file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* ============================================= */
  /*   Create and save the data attributes         */ 
  /* ============================================= */
  int attrblistLEN=13;
  // H5T_NATIVE_DOUBLE  => H5T_IEEE_F64LE
  // H5T_NATIVE_INT     => H5T_STD_I64LE
  struct {char  name[128];hid_t H5type; void * data;} attrblist[] ={
    {"poly K"      , H5T_NATIVE_DOUBLE, eos_k  },
    {"poly Gamma"  , H5T_NATIVE_DOUBLE, Gamma_P},
    {"axis_ratio"  , H5T_NATIVE_DOUBLE, r_ratio},
    {"re"          , H5T_NATIVE_DOUBLE, r_e},
    {"rhoc"        , H5T_NATIVE_DOUBLE, rho0_center},
    {"A diff"      , H5T_NATIVE_DOUBLE, A_diff},
    {"Omega"       , H5T_NATIVE_DOUBLE, Omega},
    {"Omega_e"     , H5T_NATIVE_DOUBLE, Omega_e},
    {"SDIV"        , H5T_NATIVE_INT, sdiv},
    {"MDIV"         ,H5T_NATIVE_INT, mdiv},
    {"rotation_type",H5T_C_S1,rotation_type},
    {"eos_type"     ,H5T_C_S1,eos_type},
    {"eos_file"     ,H5T_C_S1,eos_file}
  };
  /* ============================================= */
  /*  Read the data set for variables and save it  */ 
  /* ============================================= */
  int varlistLEN = 9;
  struct {char  name[128];double ** data;char description[128];} varlist[] ={
    {"/rho_potential",rho_potential,"Values for RNSID variable rho_potential"},
    {"/gama"         ,gama         ,"Values for RNSID variable gama         "},
    {"/alpha"        ,alpha        ,"Values for RNSID variable alpha        "},
    {"/omega"        ,omega        ,"Values for RNSID variable omega        "},
    {"/energy"       ,energy       ,"Values for RNSID variable energy       "},
    {"/pressure"     ,pressure     ,"Values for RNSID variable pressure     "},
    {"/enthalpy"     ,enthalpy     ,"Values for RNSID variable enthalpy     "},
    {"/velocity_sq"  ,velocity_sq  ,"Values for RNSID variable velocity_sq  "},
    {"/Omega_diff"   ,Omega_diff   ,"Values for RNSID variable Omega_diff   "}
  };


  /* ============================================= */
  /*   Create and save the data attributes         */ 
  /* ============================================= */
  for(varIndex = 0; varIndex< attrblistLEN;varIndex++) {
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
        attributeH5type = H5Tcopy(H5T_C_S1);
  	H5Tset_size(attributeH5type,strlen(attrblist[varIndex].data)+1);
  	dims[0]=1;
    } else {
      dims[0]=1;
      attributeH5type=attrblist[varIndex].H5type;
      
    }
    dataspace_id = H5Screate_simple(1, dims, NULL);
    attribute_id = H5Acreate (file_id, 
  			       attrblist[varIndex].name, attributeH5type,
                        dataspace_id,  H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id,attributeH5type, attrblist[varIndex].data);
    status = H5Aclose(attribute_id);
    status = H5Sclose(dataspace_id);
#ifdef RNS_SEQ_COMPILATION
    if (attrblist[varIndex].H5type == H5T_NATIVE_DOUBLE) {
      printf("Save attribute %s value is %5.4e\n",attrblist[varIndex].name,*((double *)attrblist[varIndex].data));
    } else if (attrblist[varIndex].H5type == H5T_NATIVE_INT) {
      printf("Save attribute %s value is %d\n",attrblist[varIndex].name,*((int *)attrblist[varIndex].data));
    } else if (attrblist[varIndex].H5type == H5T_C_S1) {
      printf("Save attribute %s value is %s\n",attrblist[varIndex].name,(char*)attrblist[varIndex].data);
    } else {
      printf("Save attribute %s value is not knowns\n",attrblist[varIndex].name);
    }
#endif

    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
      status= H5Tclose(attributeH5type);
    }
  }

  /* ============================================= */
  /*   Create and save the data set GRID points    */ 
  /* ============================================= */
  dims[0]= (*sdiv)+1;
  dataspace_id = H5Screate_simple(1,dims, NULL);
  dataset_id = H5Dcreate(file_id,"/s_qp", 
                         H5T_NATIVE_DOUBLE, dataspace_id, 
  			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  		      s_qp);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  dims[0]= (*mdiv)+1;
  dataspace_id = H5Screate_simple(1,dims, NULL);
  dataset_id = H5Dcreate(file_id,"/mu", 
                         H5T_NATIVE_DOUBLE, dataspace_id, 
  			 H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
  		      mu);
  status = H5Sclose(dataspace_id);
  status = H5Dclose(dataset_id);
  
  /* ============================================= */
  /*   Create and save the data set for variables  */ 
  /* ============================================= */
  attributeH5c128 = H5Tcopy(H5T_C_S1);
  H5Tset_size(attributeH5c128,128);
  
  for(varIndex = 0; varIndex< varlistLEN;varIndex++) {
    dset_data = malloc( *sdiv * *mdiv * sizeof(double) );
    dims[0]   = *sdiv; 
    dims[1]   = *mdiv; 
    var = varlist[varIndex].data;
    double minval,maxval;
    minval = 100.0;
    maxval =0.0;
 
    for (i = 0; i < *sdiv; i++) {
      for (j = 0; j < *mdiv; j++) {
	dset_data[i*(*mdiv)+j] = var[i+1][j+1];
	if (minval > dset_data[i*(*mdiv)+j])
          minval = dset_data[i*(*mdiv)+j];
	if (maxval < dset_data[i*(*mdiv)+j])
          maxval = dset_data[i*(*mdiv)+j];
      }
    }
#ifdef RNS_SEQ_COMPILATION
    printf("Save variable %s %g %g \n",varlist[varIndex].name,minval,maxval);
#endif
    dataspace_id = H5Screate_simple(2, dims, NULL);
    dataset_id = H5Dcreate(file_id, varlist[varIndex].name, 
                           H5T_NATIVE_DOUBLE, dataspace_id, 
			   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
		      dset_data);
    status = H5Sclose(dataspace_id);
    status = H5Dclose(dataset_id);
    free(dset_data);
    /* ========================================================= */
    /* ---- Now add the attributes to the Just created data sets */
    /* ========================================================= */
    dataset_id = H5Dopen(file_id,varlist[varIndex].name, H5P_DEFAULT);
    dims[0]=1;
    dataspace_id = H5Screate_simple(1, dims, NULL);
    attribute_id = H5Acreate (dataset_id,"Description",attributeH5c128,
			      dataspace_id,  H5P_DEFAULT, H5P_DEFAULT);
    status = H5Awrite(attribute_id,attributeH5c128, varlist[varIndex].description);
    status = H5Aclose(attribute_id);
    status = H5Dclose(dataset_id);


  }

  status= H5Tclose(attributeH5c128);

  /* Terminate access to the file. */
  status = H5Fclose(file_id); 
}


void hdf5_read_var(int *sdiv, int *mdiv,
		   char fname[],
		   /* ------------------- */
                   char eos_type[],
                   char eos_file[],
                   double *eos_k,
                   double *Gamma_P,
		   /* ------------------- */
                   char   rotation_type[],
                   double *A_diff,
		   /* ------------------- */
                   double *r_ratio,
                   double *rho0_center,
                   double *r_e,	
		   /* ------------------- */
                   double *s_qp,  /*SDIV+1*/
                   double *mu,    /*MDIV+1*/
		   /* ------------------- */
                   double **rho_potential, 
                   double **gama, 
                   double **alpha,
                   double **omega,
                   double **energy,
                   double **pressure,
                   double **enthalpy,
                   double **velocity_sq,
                   double *Omega,
                   double *Omega_e,
		   double **Omega_diff)
{

  hid_t       file_id;   /* file identifier */
  herr_t      status;

  hid_t       attribute_id; /* identifiers for attributes*/
  hid_t       attributeH5type,attributeH5c128;
  hid_t       datasetH5type;
  hid_t       dataset_id, dataspace_id;  /* identifiers for dsets*/
  hsize_t     dims[3];

  int         i,j;
  double      *dset_data;
  double      **var;
  int         varIndex;

  /* =========================================== */
  /* Create a new file using default properties. */
  /* =========================================== */
  /* file_id = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT); */
  file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);

  /* ============================================= */
  /*   Create and save the data attributes         */ 
  /* ============================================= */
  int attrblistLEN=12;
  struct {char  name[128];hid_t H5type; void * data;} attrblist[] ={
    {"poly K"      , H5T_NATIVE_DOUBLE, eos_k  },
    {"poly Gamma"  , H5T_NATIVE_DOUBLE, Gamma_P},
    {"axis_ratio"  , H5T_NATIVE_DOUBLE, r_ratio},
    {"re"          , H5T_NATIVE_DOUBLE, r_e},
    {"rhoc"        , H5T_NATIVE_DOUBLE, rho0_center},
    {"A diff"      , H5T_NATIVE_DOUBLE, A_diff},
    {"Omega"       , H5T_NATIVE_DOUBLE, Omega},
    {"Omega_e"     , H5T_NATIVE_DOUBLE, Omega_e},
    {"SDIV"        , H5T_NATIVE_INT, sdiv},
    {"MDIV"         ,H5T_NATIVE_INT, mdiv},
    {"rotation_type",H5T_C_S1,rotation_type},
    {"eos_type"     ,H5T_C_S1,eos_type},
    {"eos_file"     ,H5T_C_S1,eos_file}
  };
  /* ============================================= */
  /*  Read the data set for variables and save it  */ 
  /* ============================================= */
  int varlistLEN = 9;
  struct {char  name[128];double ** data;char description[128];} varlist[] ={
    {"/rho_potential",rho_potential,"Values for RNSID variable rho_potential"},
    {"/gama"         ,gama         ,"Values for RNSID variable gama         "},
    {"/alpha"        ,alpha        ,"Values for RNSID variable alpha        "},
    {"/omega"        ,omega        ,"Values for RNSID variable omega        "},
    {"/energy"       ,energy       ,"Values for RNSID variable energy       "},
    {"/pressure"     ,pressure     ,"Values for RNSID variable pressure     "},
    {"/enthalpy"     ,enthalpy     ,"Values for RNSID variable enthalpy     "},
    {"/velocity_sq"  ,velocity_sq  ,"Values for RNSID variable velocity_sq  "},
    {"/Omega_diff"   ,Omega_diff   ,"Values for RNSID variable Omega_diff   "}
  };

  for(varIndex = 0; varIndex< attrblistLEN;varIndex++) {
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
        attributeH5type = H5Tcopy(H5T_C_S1);
  	H5Tset_size(attributeH5type,strlen(attrblist[varIndex].data)+1);
  	dims[0]=1;
    } else {
      dims[0]=1;
      attributeH5type=attrblist[varIndex].H5type;
      
    }

    dataspace_id  = H5Screate_simple(1, dims, NULL);
    attribute_id  = H5Aopen(file_id,attrblist[varIndex].name, H5P_DEFAULT);
    attributeH5type = H5Aget_type(attribute_id);
    attributeH5type = H5Tget_native_type(attributeH5type,H5T_DIR_DEFAULT);

    /* -------------------------------------------------------------
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
    } else {
      if ( attributeH5type != attrblist[varIndex].H5type) {
	printf("Attribute type missmatch %d %d %d !\n",attributeH5type,attrblist[varIndex].H5type,H5T_IEEE_F64LE); 
	attributeH5type = attrblist[varIndex].H5type;
      }
    } 
    -------------------------------------------------------------- */
    //    attribute_id = H5Acreate (file_id, 
    //	  	       attrblist[varIndex].name, attributeH5type,
    //                  dataspace_id,  H5P_DEFAULT, H5P_DEFAULT);
    status = H5Aread(attribute_id,attributeH5type, attrblist[varIndex].data);
    status = H5Aclose(attribute_id);

#ifdef RNS_SEQ_COMPILATION
    if (attrblist[varIndex].H5type == H5T_NATIVE_DOUBLE)
      printf("Read attribute %s value is %5.4e\n",attrblist[varIndex].name,*((double *)attrblist[varIndex].data));
    else if (attrblist[varIndex].H5type == H5T_NATIVE_INT) 
      printf("Read attribute %s value is %d\n",attrblist[varIndex].name,*((int *)attrblist[varIndex].data));
    else if (attrblist[varIndex].H5type == H5T_C_S1) 
      printf("Read attribute %s value is %s\n",attrblist[varIndex].name,(char*)attrblist[varIndex].data);
    else 
      printf("Read attribute %s value is not knowns\n",attrblist[varIndex].name);
#endif

    /// status = H5Sclose(dataspace_id);
    if ( attrblist[varIndex].H5type == H5T_C_S1 ) {
      status= H5Tclose(attributeH5type);
    }
  }

  /* ============================================= */
  /*  Read the data set for GRID points            */
  /* ============================================= */
  dataset_id =H5Dopen(file_id,"/s_qp", H5P_DEFAULT);

  datasetH5type =H5Dget_type(dataset_id);
  datasetH5type =H5Tget_native_type(datasetH5type,H5T_DIR_DEFAULT); 
  status = H5Dread(dataset_id,datasetH5type, H5S_ALL, H5S_ALL, H5P_DEFAULT,s_qp);
  status = H5Dclose(dataset_id);

  dataset_id =H5Dopen(file_id,"/mu", H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,mu);
  status = H5Dclose(dataset_id);

  /* attributeH5c128 = H5Tcopy(H5T_C_S1); */
  /* H5Tset_size(attributeH5c128,128);    */
  /* status= H5Tclose(attributeH5c128);   */

  /* ============================================= */
  /*  Read the data set for variables              */
  /* ============================================= */
  for(varIndex = 0; varIndex< varlistLEN;varIndex++) {
    dset_data = malloc( *sdiv * *mdiv * sizeof(double) );
    dims[0]   = *sdiv; 
    dims[1]   = *mdiv; 
    dataset_id =H5Dopen(file_id,varlist[varIndex].name, H5P_DEFAULT);
    status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dset_data);
    //status = H5Sclose(dataspace_id);
    var = varlist[varIndex].data;
    double minval,maxval;
    minval = 100.0;
    maxval =0.0;
    for (i = 0; i < *sdiv; i++) {
      for (j = 0; j < *mdiv; j++) {
	var[i+1][j+1] = dset_data[i*(*mdiv)+j] ;
 	if (minval > dset_data[i*(*mdiv)+j])
	  minval = dset_data[i*(*mdiv)+j];
 	if (maxval < dset_data[i*(*mdiv)+j])
	  maxval = dset_data[i*(*mdiv)+j];
      }
    }
    status = H5Dclose(dataset_id);
#ifdef RNS_SEQ_COMPILATION
    printf("Read variable %s %g %g \n",varlist[varIndex].name,minval,maxval);
#endif
    free(dset_data);

  }


  /* Terminate access to the file. */
  status = H5Fclose(file_id); 
}

