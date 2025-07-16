# ISDA CDS Example

# Check if ISDA CDS support is enabled
if(TARGET QuantLib AND ISDA_CDS_LIBRARY)
    add_executable(IsdaCdsExample IsdaCdsExample.cpp)
    
    target_link_libraries(IsdaCdsExample 
        QuantLib
        ${ISDA_CDS_LIBRARY}
    )
    
    target_include_directories(IsdaCdsExample PRIVATE 
        ${ISDA_CDS_INCLUDE_DIR}
    )
    
    # Add to examples
    set_target_properties(IsdaCdsExample PROPERTIES
        FOLDER "Examples"
    )
    
    message(STATUS "IsdaCdsExample will be built")
else()
    message(STATUS "IsdaCdsExample will be skipped (ISDA CDS library not found)")
endif()