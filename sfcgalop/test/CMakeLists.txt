add_executable(test_sfcgalop
    test_minimal.cpp
    ../io.cpp
    ../util.cpp
    ../error_handler.cpp
    ../text_ui.cpp
    ../constructors.cpp
    ../operations/operations.cpp
    ../operations/operations_registry.cpp
    ../operations/operations_analysis.cpp
    ../operations/operations_collections.cpp
    ../operations/operations_construction.cpp
    ../operations/operations_constructors.cpp
    ../operations/operations_conversions.cpp
    ../operations/operations_metrics.cpp
    ../operations/operations_predicates.cpp
    ../operations/operations_set.cpp
    ../operations/operations_transformations.cpp
)

target_include_directories(test_sfcgalop PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/..
    ${SFCGAL_INCLUDE_DIRS}
)

target_link_libraries(test_sfcgalop PRIVATE
    SFCGAL
    Boost::unit_test_framework
    ${CGAL_3RD_PARTY_LIBRARIES}
)

add_test(NAME sfcgalop_unit_tests COMMAND test_sfcgalop)
