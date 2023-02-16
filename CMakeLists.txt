project(CustomThreeDeeSDL)
cmake_minimum_required(VERSION 3.23)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXE_LINKER_FLAGS "-static")

include_directories(C:/libs/SDL2/include) # C:/deps/SDL2_image/include
link_directories(C:/libs/SDL2/lib/x64) # C:/deps/SDL2_image/lib/x64

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME})
option(BUILD_CONSOLE "Build with console window" ON)

if (WIN32)
    if (BUILD_CONSOLE)
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mconsole")
    else()
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -mwindows")
    endif()
endif()

add_executable(ThreeDee src/main.cpp src/Clock.cpp src/Clock.h deps/icon_resource.rc src/Colors.cpp src/Colors.h)
target_link_libraries(ThreeDee SDL2main SDL2) # SDL2_image

function(copy_files SOURCE_DIR DEST_DIR)
    file(GLOB files "${SOURCE_DIR}/*")
    foreach(file ${files})
        get_filename_component(filename ${file} NAME)
        configure_file(${file} ${DEST_DIR}/${filename} COPYONLY)
    endforeach()
endfunction()

copy_files(${CMAKE_CURRENT_SOURCE_DIR}/bin ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME})
