#ifndef CS247_PROG_H
#define CS247_PROG_H

#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <string>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// framework includes
#include "glslprogram.h"
#include "vboquad.h"


////////////////
// Structures //
////////////////

// window size
const unsigned int gWindowWidth = 512;
const unsigned int gWindowHeight = 512;

int current_scalar_field;
int data_size;
bool en_arrow;
bool en_streamline;
bool en_pathline;

int sampling_rate;
float dt;


//////////////////////
//  Global defines  //
//////////////////////
#define TIMER_FREQUENCY_MILLIS  50

// integration methods
#define INTEG_EULER 0
#define INTEG_RK2   1
#define INTEG_RK4   2

// colormap modes
#define CMAP_OFF      0
#define CMAP_RAINBOW  1
#define CMAP_COOLWARM 2
#define CMAP_SOLID   -1   // overlay: ignore texture, use solid color

//////////////////////
// Global variables //
//////////////////////

// Handle of the window we're rendering to
static GLFWwindow* window;

char bmModifiers;	// keyboard modifiers (e.g. ctrl,...)

int clearColor;

// data handling
char* filenames[ 3 ];
bool grid_data_loaded;
bool scalar_data_loaded;
unsigned short vol_dim[ 3 ];
float* vector_array;
float* scalar_fields;
float* scalar_bounds;

GLuint scalar_field_texture;

int num_scalar_fields;
int num_timesteps; //stores number of time steps

int loaded_file;
int loaded_timestep;
float timestep;

int view_width, view_height; // height and width of entire view

GLuint displayList_idx;

int toggle_xy;

// integration / interaction state
int integration_method;     // INTEG_EULER / INTEG_RK2 / INTEG_RK4
bool en_constant_length;    // true: constant arrow length, false: speed-scaled
int rake_mode;              // 0: off, 1: vertical rake, 2: horizontal rake

#define RAKE_OFF        0
#define RAKE_VERTICAL   1
#define RAKE_HORIZONTAL 2

// colormap state
int colormap_mode;          // CMAP_OFF / CMAP_RAINBOW / CMAP_COOLWARM
float blend_factor;         // [0,1]: blend between grayscale and colormap

// stored seeds (so streamlines can be recomputed when the timestep changes)
struct Seed { float x; float y; int t; };
std::vector<Seed> streamline_seeds;
std::vector<Seed> pathline_seeds;

// flat vertex buffers (pairs of NDC coords, drawn as GL_LINES)
std::vector<float> glyph_vertices;
std::vector<float> streamline_vertices;
std::vector<float> pathline_vertices;

// GL handles for overlays
GLuint glyphVAO, glyphVBO;
GLuint streamVAO, streamVBO;
GLuint pathVAO, pathVBO;

////////////////
// Prototypes //
////////////////

void drawGlyphs();

void computeStreamline(int x, int y);

void computePathline(int x, int y, int t);

void loadNextTimestep( void );

void LoadData( char* base_filename );
void LoadVectorData( const char* filename );

void DownloadScalarFieldAsTexture( void );
void initGL( void );

void reset_rendering_props( void );

void recomputeAllStreamlines( void );
void uploadOverlay( GLuint vao, GLuint vbo, const std::vector<float>& data );
void drawOverlay( GLuint vao, GLsizei vertex_count, const glm::vec4& color );

// make quad to load texture to
VBOQuad quad;

// GLSL
GLSLProgram vectorProgram;
glm::mat4 model;


#endif //CS247_PROG_H
