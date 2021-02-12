# for Ubuntu Linux 18.04, you should perform:
# $ sudo apt install libgdal-dev gdal-bin g++
#
# --- Compilers ---
CXX  := g++

# --- Flags for compiler and linker. ---
AGAP_ENV := false
ifeq ($(AGAP_ENV), true)
  INCS := -I/usr/local/cuda/include -I/dats4/lulc/opt/saclass-ext-lib/include
  LDFLAGS := -L/dats4/lulc/opt/saclass-ext-lib/lib -Xlinker -rpath -Xlinker /dats4/lulc/opt/saclass-ext-lib/lib
else
  INCS := -I/usr/local/include -I/usr/local/cuda/include -I/usr/include/gdal
  LDFLAGS := -L/usr/local/lib
endif

CXXFLAGS := -std=c++0x -Wall

CXX_DEBUG_FLAGS := -g -O0 -ggdb3 -lefence
CXX_RELEASE_FLAGS := -O3 -fopenmp

# --- Target rules ---
all: sun dem2shadow shadow_slope_correct_Landsat5

dem2shadow:
	$(CXX) src/$@.cpp src/vector_math.cpp $(CXXFLAGS) $(CXX_RELEASE_FLAGS) $(INCS) $(LDFLAGS) $(OBJS) -lgdal -o $@

shadow_slope_correct_Landsat5:
	$(CXX) src/$@.cpp $(CXXFLAGS) $(CXX_RELEASE_FLAGS) $(INCS) $(LDFLAGS) $(OBJS) -lgdal -o $@

sun:
	$(CXX) src/$@.c $(CXXFLAGS) $(CXX_RELEASE_FLAGS) $(INCS) $(LDFLAGS) $(OBJS) -lgdal -o $@
