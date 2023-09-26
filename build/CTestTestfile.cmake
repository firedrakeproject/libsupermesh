# CMake generated Testfile for 
# Source directory: /Users/ddolci/firedrake/src/libsupermesh
# Build directory: /Users/ddolci/firedrake/src/libsupermesh/build
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_element_intersector "test_element_intersector")
set_tests_properties(test_element_intersector PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_intersection_finder_2d "test_intersection_finder_2d")
set_tests_properties(test_intersection_finder_2d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_intersection_finder_3d "test_intersection_finder_3d")
set_tests_properties(test_intersection_finder_3d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_intersection_finder_completeness_1d "test_intersection_finder_completeness_1d")
set_tests_properties(test_intersection_finder_completeness_1d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_intersection_finder_completeness_2d "test_intersection_finder_completeness_2d")
set_tests_properties(test_intersection_finder_completeness_2d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_intersection_finder_completeness_3d "test_intersection_finder_completeness_3d")
set_tests_properties(test_intersection_finder_completeness_3d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_interval_intersector "test_interval_intersector")
set_tests_properties(test_interval_intersector PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_parallel_p1_inner_product_2d "/opt/homebrew/bin/mpiexec" "-n" "4" "test_parallel_p1_inner_product_2d")
set_tests_properties(test_parallel_p1_inner_product_2d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;214;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_parallel_p1_inner_product_3d "/opt/homebrew/bin/mpiexec" "-n" "4" "test_parallel_p1_inner_product_3d")
set_tests_properties(test_parallel_p1_inner_product_3d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;214;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_parallel_p2_inner_product_2d "/opt/homebrew/bin/mpiexec" "-n" "4" "test_parallel_p2_inner_product_2d")
set_tests_properties(test_parallel_p2_inner_product_2d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;214;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_parallel_supermesh_2d "/opt/homebrew/bin/mpiexec" "-n" "4" "test_parallel_supermesh_2d")
set_tests_properties(test_parallel_supermesh_2d PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;214;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_tet_intersector "test_tet_intersector")
set_tests_properties(test_tet_intersector PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
add_test(test_tri_intersector "test_tri_intersector")
set_tests_properties(test_tri_intersector PROPERTIES  FAIL_REGULAR_EXPRESSION "Fail:" _BACKTRACE_TRIPLES "/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;216;add_test;/Users/ddolci/firedrake/src/libsupermesh/CMakeLists.txt;0;")
subdirs("spatialindex-1.8.5")
