/* 
 * File:   dem2shadow.cpp
 * Author: Kenlo Nishida Nasahara
 * mapping shadow area from DEM (digital evevation model)
 *   shadow --> given value by --valshadow= option (default=255)
 *   non-shadow --> 100*cos(angle between slope normal and sun direction)
 * Created on 2017/11/15
 * compile (on jsbach03):      g++ dem2shadow.cpp vector_math.cpp -lgdal -std=c++0x -fopenmp -o dem2shadow
 *         (on Ubuntu 16.04):  g++ dem2shadow.cpp vector_math.cpp -lgdal -std=c++11 -fopenmp -o dem2shadow
 * example: $ ./dem2shadow 67.9 60.9 N036E140_AVE_DSM.tif shadow.tif
 * example: $ ./dem2shadow 67.9 60.9 GSI_DEM10m_N36E140.tif shadow.tif --zscale=0.1 
 * note: files are in GeoTIFF (both input and output)
 * note: solar_azimuth 0 is east. 90 is south. 180 is west.
 * note: solar_zenith 0 is top. 90 is horizon.
 * note: zscale is unit of height. Default is 1. But for example, if 0.1 m unit (such as GSI 10m DSM), zscale is 0.1.
 */

#include <gdal_priv.h>
#include <cpl_conv.h> // for CPLMalloc()
#include <ogr_spatialref.h>
#include <stdint.h>
#include <cstdint>
#include <string>
#include <iostream>
#include <math.h>
#include <unistd.h>
#include <omp.h>
#include "cmdline.h"
#include "vector_math.h"
# define _USE_MATH_DEFINES      // for use of M_PI
# define EarthRadius 6371000;   // in meters
# define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
# define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
# define RETURN_ERROR 1


int main(int argc, char* argv[])
{
    cmdline::parser p;
    p.set_program_name(std::string(argv[0])+" solar_azimuth solar_zenith input(DEM)filename output_filename");
    p.add<float>("zscale", 0, "z scale of DEM. For GSI DSM, 0.1", false, 1.0);
    p.add<int>("outcol", 0, "columns (pixels) in output", false, 0);
    p.add<int>("outrow", 0, "rows (lines) in output", false, 0);
    p.add<int>("valshadow", 0, "value for a shadow pixel", false, 255);
    p.add("quiet", 'q', "display no messages");
    p.add("help", 'h', "print help");
    if (!p.parse(argc, argv)||p.exist("help")){
        std::cout<<p.error_full()<<p.usage();
        return RETURN_ERROR;
    }
    
    double dzDEM=p.get<float>("zscale");  // lattice space in z direction of original image (DEM); normally 1.0 m
    int32_t nXSizeOUT=p.get<int>("outcol");
    int32_t nYSizeOUT=p.get<int>("outrow");
    int valshadow=p.get<int>("valshadow");
    
    if (argc<4){
        printf("Error in %s: Incomplete arguments! Check $ %s --help\n", argv[0], argv[0]);
        return RETURN_ERROR;
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  I got solar azimuth=%s degree; solar zenith=%s degree; DEM: %s; output: %s\n", argv[0], argv[1], argv[2], argv[3], argv[4]);

    double solaz=atof(argv[1])*M_PI/180;   // solar azimuth angle
    double solzn=atof(argv[2])*M_PI/180;   // solar zenith angle

// initialize GDAL interface    
    GDALAllRegister();     // register GDAL drivers
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");

// Read Digital Elevation Model (DEM) file
    GDALDataset  *poSrcDS;    // pointer to source dataset
    if ((poSrcDS = (GDALDataset *) GDALOpen(argv[3], GA_ReadOnly ))==NULL) exit(1);  // Open the GeoTIFF file. If error, exit!
   // Get spatial reference system
    const char *pszWkt = poSrcDS->GetProjectionRef();
    OGRSpatialReference oSRS;
    oSRS.importFromWkt((char **)&pszWkt);
    int NumBand=poSrcDS->GetRasterCount();
   // Get affine transformation parameters
    double        adfGeoTransform[6];   // array of double float for affine transformation parameters
    if( poSrcDS->GetGeoTransform( adfGeoTransform ) == CE_None ){
    }
  // Get single-band raster data body
    GDALRasterBand  *poBand;
    int             nBlockXSize, nBlockYSize;
    int             bGotMin, bGotMax;
    double          adfMinMax[2];
    poBand = poSrcDS->GetRasterBand( 1 );
    poBand->GetBlockSize( &nBlockXSize, &nBlockYSize );
    adfMinMax[0] = poBand->GetMinimum( &bGotMin );
    adfMinMax[1] = poBand->GetMaximum( &bGotMax );
    if( ! (bGotMin && bGotMax) )
    GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
    if( poBand->GetOverviewCount() > 0 ) printf( "Band has %d overviews.\n", poBand->GetOverviewCount() );
    if( poBand->GetColorTable() != NULL ) printf( "Band has a color table with %d entries.\n", poBand->GetColorTable()->GetColorEntryCount() );

    int32_t   nXSizeDEM = poBand->GetXSize();
    int32_t   nYSizeDEM = poBand->GetYSize();
    if (nXSizeOUT<1) nXSizeOUT=nXSizeDEM;
    if (nYSizeOUT<1) nYSizeOUT=nYSizeDEM;
    
    GInt16 *dem;    // array to store DEM data
    if ((dem = (GInt16 *) CPLMalloc(sizeof(GInt16)*nXSizeDEM*nYSizeDEM))==NULL){    // memory allocation. If fail, give error message and exit.
        printf("Error in %s: Failed to allocate memory for DEM!\n", argv[0]);
        return RETURN_ERROR;
    }
    if(poBand->RasterIO( GF_Read, 0, 0, nXSizeDEM, nYSizeDEM,dem, nXSizeDEM, nYSizeDEM, GDT_Int16,0, 0 )==3){
        printf("Error in %s: Failed reading DEM!\n", argv[0]);
        return RETURN_ERROR;
    }   // read DEM data
    if (!p.exist("quiet"))
      printf("[Runnning %s]  Reading %s successful!\n", argv[0], argv[3]);

  // adjustment of anisotropy of pixel lattice
    double aboutLatitude=adfGeoTransform[3]*M_PI/180;    // approximate latitude to estimate the lattice space in x (west-east) direction.
    double dxDEM=fabs(adfGeoTransform[1])*(M_PI/180)*cos(aboutLatitude)*EarthRadius;  // lattice space in x direction of original image (DEM)
    double dyDEM=fabs(adfGeoTransform[5])*(M_PI/180)*EarthRadius;                     // lattice space in y direction of original image (DEM)
    if (dxDEM == 0 || dyDEM == 0 || dzDEM == 0){
        printf("[Error in %s]  Lattice space (%lf, %lf, %lf) is strange!", argv[0], dxDEM, dyDEM, dzDEM);
        return RETURN_ERROR;
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  DEM pixels: (%d, %d);  dxDEM=%.3f, dyDEM=%.3f, dzDEM=%.3f; output pixels: (%d, %d); shadow value=%d\n", 
                  argv[0], nXSizeDEM, nYSizeDEM, dxDEM, dyDEM, dzDEM, nXSizeOUT, nYSizeOUT, valshadow);

// image processing
    
   // create and initialize an array "shadowmap" to store shadow (1) and non-shadow (0) pixels.
    unsigned char *shadowmap;    
    if ((shadowmap = (unsigned char *) CPLMalloc(sizeof(unsigned char)*nXSizeOUT*nYSizeOUT))==NULL){  // memory allocation
        printf("Error in %s: Failed to allocate memory for shadow map!\n", argv[0]);
        return RETURN_ERROR;
    }

    VECT sunvectRW=vectdef(sin(solzn)*cos(solaz), sin(solzn)*sin(solaz), cos(solzn));   // vector to the sun in usual (real-world; isotropic) coordinate 
    int32_t ox, oy; // index for output image
   // calculating the gradient vector and slope at each pixel 
    for (oy=0; oy<nYSizeOUT; oy++){
       for (ox=0; ox<nXSizeOUT; ox++){
           double oxc=ox+0.5;   // "c" means center
           double oyc=oy+0.5;   // "c" means center
           int32_t mxc=( oxc   *nXSizeDEM/nXSizeOUT);  // "m" means DEM; "c" means center
           int32_t myc=( oyc   *nYSizeDEM/nYSizeOUT);
           int32_t mxw=((oxc-1)*nXSizeDEM/nXSizeOUT); if (mxw<0) mxw=0;                     // "w" means west
           int32_t myn=((oyc-1)*nYSizeDEM/nYSizeOUT); if (myn<0) myn=0;                     // "n" means north
           int32_t mxe=((oxc+1)*nXSizeDEM/nXSizeOUT); if (nXSizeDEM<=mxe) mxw=nXSizeDEM-1;  // "e" means east
           int32_t mys=((oyc+1)*nYSizeDEM/nYSizeOUT); if (nYSizeDEM<=mys) mys=nYSizeDEM-1;  // "s" means south
           double dzpdx= dzDEM * (dem[mxe+myc*nXSizeDEM] - dem[mxw+myc*nXSizeDEM]) / (dxDEM*(mxe-mxw));   // dz/dx
           double dzpdy= dzDEM * (dem[mxc+mys*nXSizeDEM] - dem[mxc+myn*nXSizeDEM]) / (dyDEM*(mys-myn));   // dzx/dy
           VECT slpvectRW=vectdef(-1*dzpdx, -1*dzpdy, 1);
           slpvectRW=slpvectRW/norm(slpvectRW);       
           shadowmap[ox+oy*nXSizeOUT]=(int)(scprod(slpvectRW, sunvectRW)*100.0);   // initialize
       }
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  gradient calculation completed.\n", argv[0]);
    
   // calculating direction of the ray tracing (sun direction)
    VECT sunvectLT=vectdef(sin(solzn)*cos(solaz)/dxDEM, sin(solzn)*sin(solaz)/dyDEM, cos(solzn)/dzDEM);  // vector to the sun in lattice coordinate (maybe anisotropic)
    sunvectLT=sunvectLT/sqrt(pow(sunvectLT.x, 2.0)+pow(sunvectLT.y, 2.0));   // horizontal length of the sun vector should be equal to one pixel.

   // setting boundaries containing the landscape (the end of ray tracing)
    int32_t mxmin=0, mxmax=nXSizeDEM, mymin=0, mymax=nYSizeDEM;
    double  mzmin=0, mzmax=0;
    int32_t mx, my; // index for DEM; "m" means DEM
    for (my=1; my<nYSizeDEM-1; my++)
       for (mx=1; mx<nXSizeDEM-1; mx++) 
           if(mzmax<dem[mx+my*nXSizeDEM]) mzmax=dem[mx+my*nXSizeDEM];  // the roof

   // scan all the pixels on the DEM and check whether each pixel is in the shadow or not. (ox, oy) means "original" x and y.
    for (oy=1; oy<nYSizeOUT-1; oy++){
      // processing message (not essential)
       if (!p.exist("quiet"))
         if (oy%(nYSizeOUT/10)==0) printf("[Runnning %s]  Calculation %d percent completed!\n", argv[0], 100*oy/nYSizeOUT);
       
       #pragma omp parallel
       {
       #pragma omp for
       for (ox=1; ox<nXSizeOUT-1; ox++){
           int32_t mx= (((double)ox + 0.5) * nXSizeDEM / nXSizeOUT);  // "m" means DEM;
           int32_t my= (((double)oy + 0.5) * nYSizeDEM / nYSizeOUT);
           double  mz=dem[mx + my * nXSizeDEM];
           VECT morig=vectdef(mx, my, mz);    // original point in DEM
           
         // estimate length of a single ray tracing
           double t_max_candidate[5] = {0,0,0,0,0};
           t_max_candidate[0]=(mxmin-mx)/sunvectLT.x;
           t_max_candidate[1]=(mxmax-mx)/sunvectLT.x;
           t_max_candidate[2]=(mymin-my)/sunvectLT.y;
           t_max_candidate[3]=(mymax-my)/sunvectLT.y;
           t_max_candidate[4]=(mzmax-mz)/sunvectLT.z;
           double t_max=10000000000000000000000.0;
           int i=0;
           for (i=0; i<5; i++) 
               if (0<=t_max_candidate[i] && t_max_candidate[i]<t_max)
                   t_max=t_max_candidate[i];   // length of ray tracing is decided by the nearest boundary in sun direction.
           
         // ray tracing
           int t;  
           for (t=0; t<t_max; t++){
               VECT rayvect=morig+sunvectLT*t;        // current position of ray along the ray trace, in the DEM 3D lattice
               int32_t mxray=(int32_t)rayvect.x;   // rounded position (suffix) of ray, in the DEM 3D lattice
               int32_t myray=(int32_t)rayvect.y;
               if (1<mxray && mxray<nXSizeDEM-1 && 1<myray && myray<nYSizeDEM-1){
                   double dzpdx=(dem[(mxray+1)+myray*nXSizeDEM]-dem[(mxray-1)+myray*nXSizeDEM])/2.0;      // partial derivative dz/dx
                   double dzpdy=(dem[mxray+(myray+1)*nXSizeDEM]-dem[mxray+(myray-1)*nXSizeDEM])/2.0;      // partial derivative dz/dy
                   double dem_height=dem[mxray+myray*nXSizeDEM]+dzpdx*(rayvect.x-mxray)+dzpdy*(rayvect.y-myray);  // DEM height above or below the ray, interpolated.
                   if (rayvect.z < dem_height){
                     shadowmap[ox+oy*nXSizeOUT]=valshadow;    // (ox, oy) is in the shadow!    
                     t=t_max+1;    // escape from the loop
                  }
               }
           }
        }
       }
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  All calculation completed!\n", argv[0]);
        
// Write data
    
GDALDataset *poDstDS;   // pointer to destination dataset
char **papszOptions = NULL;  // pointer of array of pointer to a string terminated by zero
poDstDS = poDriver->Create( argv[4], nXSizeOUT, nYSizeOUT, 1, GDT_Byte, papszOptions );
adfGeoTransform[1]=adfGeoTransform[1]*nXSizeDEM/nXSizeOUT;
adfGeoTransform[5]=adfGeoTransform[5]*nYSizeDEM/nYSizeOUT;

poDstDS->SetGeoTransform( adfGeoTransform );
oSRS.exportToWkt((char **)&pszWkt );
poDstDS->SetProjection( pszWkt );
char *pszSRS_WKT = NULL;
GByte abyRaster[512*512];
CPLFree( pszSRS_WKT );
poBand = poDstDS->GetRasterBand(1);
if ((poBand->RasterIO( GF_Write, 0, 0, nXSizeOUT, nYSizeOUT, shadowmap, nXSizeOUT, nYSizeOUT, GDT_Byte, 0, 0 ))==3){
        printf("Error in %s: Failed writing result!\n", argv[0]);
        return RETURN_ERROR;
    }
CPLFree(dem);
CPLFree(shadowmap);
GDALClose(poSrcDS);
GDALClose(poDstDS);
}
