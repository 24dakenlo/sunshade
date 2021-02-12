/* 
 * File:   dem2slope.cpp
 * Author: Kenlo Nishida Nasahara
 * mapping slope from DEM (digital evevation model)
 * Created on 2017/11/15
 * compile (ubuntu16.04): g++ dem2slope.cpp -lgdal -std=c++11 -fopenmp -o dem2slope
 * compile (jsbach03):    g++ dem2slope.cpp -lgdal -std=c++0x -fopenmp -o dem2slope
 * example: $ ./dem2slope GSI_DEM10m_N36E140.tif slope.tif --zscale=0.1
 * example: $ ./dem2slope N036E140_AVE_DSM.tif slope.tif
 * note: files are in GeoTIFF (both input and output)
 * note: zscale is unit of height. Normally it is 1. But for example, if 0.1 m unit (such as GSI 10m DSM), zscale is 0.1.
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
# define _USE_MATH_DEFINES      // for use of M_PI
# define EarthRadius 6371000;   // in meters
# define max(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b ? _a : _b; })
# define min(a,b) ({ __typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a < _b ? _a : _b; })
# define RETURN_ERROR 1


int main(int argc, char* argv[])
{
    cmdline::parser p;
    p.set_program_name(std::string(argv[0])+" input(DEM)filename output_filename");
    p.add<float>("zscale", 0, "z scale of DEM. For GSI DSM, 0.1", false, 1.0);
    p.add<int>("outcol", 0, "columns (pixels) in output", false, 0);
    p.add<int>("outrow", 0, "rows (lines) in output", false, 0);
    p.add("cos", 'c', "described in cosine * 100 (default is in degree)");
    p.add("quiet", 'q', "display no messages");
    p.add("help", 'h', "print help");
    if (!p.parse(argc, argv)||p.exist("help")){
        std::cout<<p.error_full()<<p.usage();
        return RETURN_ERROR;
    }
    
    double dzDEM=p.get<float>("zscale");  // lattice space in z direction of original image (DEM); normally 1.0 m
    int32_t nXSizeOUT=p.get<int>("outcol");
    int32_t nYSizeOUT=p.get<int>("outrow");
    
    if (argc<2){
        printf("Error in %s: Incomplete arguments! Check $ %s --help\n", argv[0], argv[0]);
        return RETURN_ERROR;
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  I got DEM: %s; output: %s\n", argv[0], argv[1], argv[2]);

// initialize GDAL interface    
    GDALAllRegister();     // register GDAL drivers
    GDALDriver *poDriver;
    poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");

// Read Digital Elevation Model (DEM) file
    GDALDataset  *poSrcDS;    // pointer to source dataset
    if ((poSrcDS = (GDALDataset *) GDALOpen(argv[1], GA_ReadOnly ))==NULL) exit(1);  // Open the GeoTIFF file. If error, exit!
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
      printf("[Runnning %s]  DEM pixels: (%d, %d);  dxDEM=%.3f, dyDEM=%.3f, dzDEM=%.3f; output pixels: (%d, %d)\n", 
                  argv[0], nXSizeDEM, nYSizeDEM, dxDEM, dyDEM, dzDEM, nXSizeOUT, nYSizeOUT);
    if (!p.exist("quiet"))
       if (p.exist("cos")) 
           printf("[Runnning %s]  slope is described in cosine * 100 in one byte\n", argv[0]);
       else 
           printf("[Runnning %s]  slope is described in degree in one byte\n", argv[0]);
    
// image processing
    
   // create and initialize an array "slopemap" to store slope pixels.
    unsigned char *slopemap;
    if ((slopemap = (unsigned char *) CPLMalloc(sizeof(unsigned char)*nXSizeOUT*nYSizeOUT))==NULL){  // memory allocation
        printf("Error in %s: Failed to allocate memory for the slope map!\n", argv[0]);
        return RETURN_ERROR;
    }

    int32_t ox, oy; // index for output image
   // calculating the gradient vector and slope at each pixel 
    for (oy=0; oy<nYSizeOUT; oy++){
       #pragma omp parallel
       {
       #pragma omp for
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
           double cos = 1/((dzpdx * dzpdx) + (dzpdy * dzpdy) + 1);
           int ans=cos*100;
           if (!p.exist("cos"))
               ans=acos(cos)*180/M_PI;
           slopemap[ox+oy*nXSizeOUT]=ans;
       }
       }
    }
    if (!p.exist("quiet"))
      printf("[Runnning %s]  All calculation completed!\n", argv[0]);
// Write data
    
GDALDataset *poDstDS;   // pointer to destination dataset
char **papszOptions = NULL;  // pointer of array of pointer to a string terminated by zero
//poDstDS = poDriver->CreateCopy( "a.tif", poSrcDS, FALSE,NULL, NULL, dem);
poDstDS = poDriver->Create( argv[2], nXSizeOUT, nYSizeOUT, 1, GDT_Byte, papszOptions );

adfGeoTransform[1]=adfGeoTransform[1]*nXSizeDEM/nXSizeOUT;
adfGeoTransform[5]=adfGeoTransform[5]*nYSizeDEM/nYSizeOUT;
poDstDS->SetGeoTransform( adfGeoTransform );
oSRS.exportToWkt((char **)&pszWkt );
poDstDS->SetProjection( pszWkt );
char *pszSRS_WKT = NULL;
GByte abyRaster[512*512];
CPLFree( pszSRS_WKT );
poBand = poDstDS->GetRasterBand(1);
if ((poBand->RasterIO( GF_Write, 0, 0, nXSizeOUT, nYSizeOUT, slopemap, nXSizeOUT, nYSizeOUT, GDT_Byte, 0, 0 ))==3){
        printf("Error in %s: Failed writing result!\n", argv[0]);
        return RETURN_ERROR;
    }
CPLFree(dem);
CPLFree(slopemap);
GDALClose(poSrcDS);
GDALClose(poDstDS);
}
