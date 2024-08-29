add_executable(BSplineCurve BSplineCurve.cpp)
target_link_libraries(BSplineCurve ql_library ${QL_THREAD_LIBRARIES})
if (QL_INSTALL_EXAMPLES)
    install(TARGETS BSplineCurve RUNTIME DESTINATION ${QL_INSTALL_EXAMPLESDIR})
endif()

