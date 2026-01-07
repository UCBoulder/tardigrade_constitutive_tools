if(${TARDIGRADE_HEADER_ONLY})
    add_library(${PROJECT_NAME} INTERFACE)
    target_include_directories(
        ${PROJECT_NAME}
        INTERFACE
            $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/${CPP_SRC_PATH}>
            $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )

    install(
        TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

    foreach(package ${ADDITIONAL_HEADER_ONLY_LIBRARIES})
        install(
            TARGETS ${package}
            EXPORT ${package}_Targets
            INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
            ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
        install(FILES ${package}.h ${package}.cpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
    endforeach(package)
else()
    add_library(${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h" "${PROJECT_NAME}.tpp")
    set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER ${PROJECT_NAME}.h)
    target_link_libraries(${PROJECT_NAME} tardigrade_error_tools)
    target_compile_options(${PROJECT_NAME} PUBLIC)

    # Local builds of upstream projects require local include paths
    if(NOT cmake_build_type_lower STREQUAL "release")
        target_include_directories(
            ${PROJECT_NAME}
            PUBLIC
                ${tardigrade_error_tools_SOURCE_DIR}/${CPP_SRC_PATH}
                ${tardigrade_vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}
        )
    endif()

    install(
        TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}_Targets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )
endif()

install(FILES ${PROJECT_NAME}.cpp ${PROJECT_NAME}.tpp DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

get_target_property(target_sources ${PROJECT_NAME} SOURCES)
foreach(source ${target_sources})
    set_property(GLOBAL APPEND PROPERTY CLANG_FORMAT_SOURCE_FILES "${CMAKE_CURRENT_SOURCE_DIR}/${source}")
endforeach(source)
