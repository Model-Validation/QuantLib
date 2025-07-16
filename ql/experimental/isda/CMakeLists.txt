# ISDA CDS Interface module for QuantLib

# Add ISDA CDS library path
set(ISDA_CDS_ROOT "C:/ws/libraries/isda_cds_model_c_v1.8.3" CACHE PATH "Path to ISDA CDS library")

# Find ISDA CDS library
find_library(ISDA_CDS_LIBRARY
    NAMES cds
    PATHS ${ISDA_CDS_ROOT}/build/Release
    NO_DEFAULT_PATH
)

# Set include directories
set(ISDA_CDS_INCLUDE_DIR ${ISDA_CDS_ROOT}/lib/include)

if(ISDA_CDS_LIBRARY AND EXISTS ${ISDA_CDS_INCLUDE_DIR})
    message(STATUS "Found ISDA CDS library: ${ISDA_CDS_LIBRARY}")
    message(STATUS "ISDA CDS include dir: ${ISDA_CDS_INCLUDE_DIR}")
    
    # Add ISDA CDS interface sources
    set(ISDA_SOURCES
        isdacdsinterface.cpp
    )
    
    # Add to parent's source list
    set(ISDA_SOURCES ${ISDA_SOURCES} PARENT_SCOPE)
    
    # Add include directories for ISDA CDS
    target_include_directories(${PROJECT_NAME} PRIVATE ${ISDA_CDS_INCLUDE_DIR})
    
    # Link ISDA CDS library
    target_link_libraries(${PROJECT_NAME} PRIVATE ${ISDA_CDS_LIBRARY})
    
    # Add preprocessor definition to indicate ISDA support
    target_compile_definitions(${PROJECT_NAME} PRIVATE QL_ENABLE_ISDA_CDS)
    
else()
    message(STATUS "ISDA CDS library not found. ISDA interface will be disabled.")
    message(STATUS "  Looked for library in: ${ISDA_CDS_ROOT}/build/Release")
    message(STATUS "  Looked for headers in: ${ISDA_CDS_INCLUDE_DIR}")
endif()