set(PYTHON_SOURCE_FILES
    pyproject.toml.in
    ${PROJECT_NAME}/_settings.py.in
    ${PROJECT_NAME}/__init__.py
    ${PROJECT_NAME}/conftest.py
    ${PROJECT_NAME}/tests/test_${PROJECT_NAME}.py
)
configure_file(${PROJECT_NAME}/_settings.py.in ${PROJECT_NAME}/_settings.py)
configure_file(pyproject.toml.in pyproject.toml)

add_custom_target(${PYTHON_TARGET} ALL
                  DEPENDS ${CYTHON_TARGET})
add_custom_command(TARGET ${PYTHON_TARGET} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   "${PROJECT_SOURCE_DIR}/${PYTHON_SRC_PATH}" "${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}"
                   COMMAND ${Python_EXECUTABLE} -m build --no-isolation
                   DEPENDS ${PYTHON_SOURCE_FILES} ${CYTHON_TARGET}
                   COMMENT "Building the ${PROJECT_NAME} Cython extension - python target"
                  )
# Add pytests as a ctest function for automated testing under unified CMake/CTest tools
if(${APPLE})
    set(LIBRARY_MODIFICATION "DYLD_FALLBACK_LIBRARY_PATH=${CMAKE_BINARY_DIR}/${CPP_SRC_PATH}")
else()
    set(LIBRARY_MODIFICATION "LD_LIBRARY_PATH=ENV{LD_LIBRARY_PATH}:${CMAKE_BINARY_DIR}/${CPP_SRC_PATH}")
endif()
add_test(
    NAME pytest
    COMMAND ${CMAKE_COMMAND} -E env ${LIBRARY_MODIFICATION} pytest
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}
)
