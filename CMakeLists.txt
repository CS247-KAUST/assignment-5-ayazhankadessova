cmake_minimum_required (VERSION 3.21)
project (assignment_5B)

# --- 1. C++ STANDARD CONFIGURATION ---
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --- 2. AUTOMATIC DEPENDENCY DOWNLOAD (FETCHCONTENT) ---
include(FetchContent)

# --- GLFW ---
FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG        3.4
)
FetchContent_MakeAvailable(glfw)

# --- GLM (Math Library) ---
FetchContent_Declare(
    glm
    GIT_REPOSITORY https://github.com/g-truc/glm.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(glm)

# --- Dear ImGui (UI Library) ---
FetchContent_Declare(
    imgui
    GIT_REPOSITORY https://github.com/ocornut/imgui.git
    GIT_TAG        docking
)
FetchContent_MakeAvailable(imgui)

# --- 3. INCLUDE DIRECTORIES ---
include_directories(common common/glad)

# OpenGL (cross-platform)
find_package(OpenGL REQUIRED)

# Adding compile units to executable
add_executable(${PROJECT_NAME}
	src/glad.c 
	src/CS247_prog.cpp
	src/glslprogram.cpp
	src/vboquad.cpp
	# --- ImGui Source Files ---
    ${imgui_SOURCE_DIR}/imgui.cpp
    ${imgui_SOURCE_DIR}/imgui_draw.cpp
    ${imgui_SOURCE_DIR}/imgui_widgets.cpp
    ${imgui_SOURCE_DIR}/imgui_tables.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
    ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl3.cpp
)


# Set project assignment as startup project in Visual Studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT assignment_3)

# C++ Compiler options
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")

# Add imgui source dir to our include path
target_include_directories(${PROJECT_NAME} PRIVATE ${imgui_SOURCE_DIR})

# --- 5. LINK LIBRARIES ---
target_link_libraries(${PROJECT_NAME} glfw OpenGL::GL glm::glm)

# imgui backends include path (needed even though we don't use imgui in code)
target_include_directories(${PROJECT_NAME} PRIVATE ${imgui_SOURCE_DIR}/backends)