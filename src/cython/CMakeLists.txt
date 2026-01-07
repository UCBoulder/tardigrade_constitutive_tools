set(CYTHON_SOURCE_FILES
    setup.py.in
    main.pyx
    tardigrade_${PROJECT_NAME}_python.pxd
    tardigrade_${PROJECT_NAME}_python.pyx
)
configure_file(setup.py.in setup.py)

add_custom_target(${CYTHON_TARGET} ALL DEPENDS ${PROJECT_NAME})
add_custom_command(
    TARGET ${CYTHON_TARGET}
    PRE_BUILD
    COMMAND
        ${CMAKE_COMMAND} -E copy_directory "${PROJECT_SOURCE_DIR}/${CYTHON_SRC_PATH}"
        "${PROJECT_BINARY_DIR}/${CYTHON_SRC_PATH}"
    COMMAND
        python setup.py build_ext --build-lib ${PROJECT_BINARY_DIR}/${PYTHON_SRC_PATH}/${PROJECT_NAME}/ --build-temp .
        DEPENDS ${CYTHON_SOURCE_FILES}
    COMMENT "BUILDING THE ${PROJECT_NAME} Cython extension - cython target"
)

install(
    FILES
        "${PROJECT_SOURCE_DIR}/${CYTHON_SRC_PATH}/${PROJECT_NAME}_python.pyx"
        "${PROJECT_SOURCE_DIR}/${CYTHON_SRC_PATH}/${PROJECT_NAME}_python.pxd"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
