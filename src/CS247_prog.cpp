// CS 247 - Scientific Visualization, KAUST
//
// Programming Assignment #5
#include <cstring>
#include <algorithm>
#include <sys/stat.h>
#include "CS247_prog.h"

// search for an asset file relative to common build-dir depths; returns the
// first existing path. Falls back to the original input.
static std::string resolveAsset(const std::string& rel)
{
    const char* prefixes[] = {
        "", "./", "../", "../../", "../../../", "../../../../"
    };
    struct stat st;
    for (size_t i = 0; i < sizeof(prefixes) / sizeof(prefixes[0]); ++i) {
        std::string p = std::string(prefixes[i]) + rel;
        if (stat(p.c_str(), &st) == 0) return p;
    }
    return rel;
}

// ---------------------------------------------------------------------------
// Small helpers: indexing, interpolation, integration
// ---------------------------------------------------------------------------

static inline bool insideGrid(float x, float y)
{
    return x >= 0.0f && x <= (float)(vol_dim[0] - 1)
        && y >= 0.0f && y <= (float)(vol_dim[1] - 1);
}

// raw vector at grid cell (x,y) at timestep t -- zero outside the grid
static inline glm::vec2 vectorAtCell(int x, int y, int t)
{
    if (x < 0 || x >= vol_dim[0] || y < 0 || y >= vol_dim[1]) return glm::vec2(0.0f);
    if (t < 0) t = 0;
    if (t >= num_timesteps) t = num_timesteps - 1;
    int slice_offset = 3 * data_size * t;
    int idx = 3 * (y * vol_dim[0] + x);
    return glm::vec2(vector_array[slice_offset + idx + 0],
                     vector_array[slice_offset + idx + 1]);
}

// bilinear interpolation of vector at fractional grid position (x,y), timestep t
static glm::vec2 sampleVectorBilinear(float x, float y, int t)
{
    int x0 = (int)std::floor(x);
    int y0 = (int)std::floor(y);
    int x1 = x0 + 1;
    int y1 = y0 + 1;
    float tx = x - (float)x0;
    float ty = y - (float)y0;

    glm::vec2 v00 = vectorAtCell(x0, y0, t);
    glm::vec2 v10 = vectorAtCell(x1, y0, t);
    glm::vec2 v01 = vectorAtCell(x0, y1, t);
    glm::vec2 v11 = vectorAtCell(x1, y1, t);

    glm::vec2 v0 = (1.0f - tx) * v00 + tx * v10;
    glm::vec2 v1 = (1.0f - tx) * v01 + tx * v11;
    return (1.0f - ty) * v0 + ty * v1;
}

// trilinear interpolation: bilinear in space + linear in time
static glm::vec2 sampleVectorTrilinear(float x, float y, float ft)
{
    int t0 = (int)std::floor(ft);
    int t1 = t0 + 1;
    float a = ft - (float)t0;

    if (t0 < 0) { t0 = 0; t1 = 0; a = 0.0f; }
    if (t1 >= num_timesteps) { t1 = num_timesteps - 1; t0 = t1; a = 0.0f; }

    glm::vec2 v0 = sampleVectorBilinear(x, y, t0);
    glm::vec2 v1 = sampleVectorBilinear(x, y, t1);
    return (1.0f - a) * v0 + a * v1;
}

// integration step on a steady (streamline) field
static glm::vec2 streamStep(glm::vec2 p, float h, int t)
{
    if (integration_method == INTEG_EULER) {
        glm::vec2 k1 = sampleVectorBilinear(p.x, p.y, t);
        return p + h * k1;
    } else if (integration_method == INTEG_RK2) {
        glm::vec2 k1 = sampleVectorBilinear(p.x, p.y, t);
        glm::vec2 k2 = sampleVectorBilinear(p.x + h * k1.x, p.y + h * k1.y, t);
        return p + 0.5f * h * (k1 + k2);
    } else { // RK4
        glm::vec2 k1 = sampleVectorBilinear(p.x, p.y, t);
        glm::vec2 k2 = sampleVectorBilinear(p.x + 0.5f * h * k1.x, p.y + 0.5f * h * k1.y, t);
        glm::vec2 k3 = sampleVectorBilinear(p.x + 0.5f * h * k2.x, p.y + 0.5f * h * k2.y, t);
        glm::vec2 k4 = sampleVectorBilinear(p.x + h * k3.x, p.y + h * k3.y, t);
        return p + (h / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
    }
}

// integration step on a time-dependent (pathline) field; advances time too
static glm::vec2 pathStep(glm::vec2 p, float ft, float h)
{
    if (integration_method == INTEG_EULER) {
        glm::vec2 k1 = sampleVectorTrilinear(p.x, p.y, ft);
        return p + h * k1;
    } else if (integration_method == INTEG_RK2) {
        glm::vec2 k1 = sampleVectorTrilinear(p.x, p.y, ft);
        glm::vec2 k2 = sampleVectorTrilinear(p.x + h * k1.x, p.y + h * k1.y, ft + h);
        return p + 0.5f * h * (k1 + k2);
    } else { // RK4
        glm::vec2 k1 = sampleVectorTrilinear(p.x,                p.y,                ft);
        glm::vec2 k2 = sampleVectorTrilinear(p.x + 0.5f * h * k1.x, p.y + 0.5f * h * k1.y, ft + 0.5f * h);
        glm::vec2 k3 = sampleVectorTrilinear(p.x + 0.5f * h * k2.x, p.y + 0.5f * h * k2.y, ft + 0.5f * h);
        glm::vec2 k4 = sampleVectorTrilinear(p.x + h * k3.x,        p.y + h * k3.y,        ft + h);
        return p + (h / 6.0f) * (k1 + 2.0f * k2 + 2.0f * k3 + k4);
    }
}

// convert grid coords to OpenGL NDC ([-1,1]^2)
static inline glm::vec2 gridToNDC(float gx, float gy)
{
    float nx = 2.0f * gx / (float)(vol_dim[0] - 1) - 1.0f;
    float ny = 2.0f * gy / (float)(vol_dim[1] - 1) - 1.0f;
    return glm::vec2(nx, ny);
}

static inline void appendSegment(std::vector<float>& out, glm::vec2 a, glm::vec2 b)
{
    out.push_back(a.x); out.push_back(a.y);
    out.push_back(b.x); out.push_back(b.y);
}

// ---------------------------------------------------------------------------
// Overlay VAO/VBO helpers
// ---------------------------------------------------------------------------

static void initOverlay(GLuint& vao, GLuint& vbo)
{
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    // 2 floats per vertex (x,y in NDC); location 0 = pos (declared vec3 in shader, z=0 implicit)
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void uploadOverlay(GLuint vao, GLuint vbo, const std::vector<float>& data)
{
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float),
                 data.empty() ? nullptr : data.data(), GL_DYNAMIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void drawOverlay(GLuint vao, GLsizei vertex_count, const glm::vec4& color)
{
    if (vertex_count <= 0) return;
    vectorProgram.setUniform("colormapMode", (int)CMAP_SOLID);
    vectorProgram.setUniform("vertexColor", color);
    vectorProgram.setUniform("model", glm::mat4(1.0f));
    glBindVertexArray(vao);
    glDrawArrays(GL_LINES, 0, vertex_count);
    glBindVertexArray(0);
}

// ---------------------------------------------------------------------------
// Streamline / pathline integration
// ---------------------------------------------------------------------------

static const int   MAX_STEPS       = 2000;
static const float MAX_LENGTH      = 400.0f;   // accumulated length in grid units
static const float MIN_VEC_MAG     = 1.0e-5f;

// integrate a streamline starting at p0, in given direction (+1 or -1), append
// segments (in NDC) to the output buffer.
static void integrateStream(glm::vec2 p0, int t, float dir,
                            std::vector<float>& out)
{
    glm::vec2 p = p0;
    float h = dt * dir;
    float length = 0.0f;
    glm::vec2 prev_ndc = gridToNDC(p.x, p.y);

    for (int i = 0; i < MAX_STEPS; ++i) {
        // stopping: zero vector
        glm::vec2 v = sampleVectorBilinear(p.x, p.y, t);
        if (glm::length(v) < MIN_VEC_MAG) break;

        glm::vec2 np = streamStep(p, h, t);

        // stopping: boundary
        if (!insideGrid(np.x, np.y)) break;

        float seg_len = glm::length(np - p);
        length += seg_len;
        if (length > MAX_LENGTH) break;

        glm::vec2 cur_ndc = gridToNDC(np.x, np.y);
        appendSegment(out, prev_ndc, cur_ndc);
        prev_ndc = cur_ndc;
        p = np;
    }
}

// integrate a pathline. forward: time advances with +dt, backward: -dt
static void integratePath(glm::vec2 p0, float t0, float dir,
                          std::vector<float>& out)
{
    glm::vec2 p = p0;
    float ft = t0;
    float h = dt * dir;
    float length = 0.0f;
    glm::vec2 prev_ndc = gridToNDC(p.x, p.y);

    for (int i = 0; i < MAX_STEPS; ++i) {
        // stopping: zero vector
        glm::vec2 v = sampleVectorTrilinear(p.x, p.y, ft);
        if (glm::length(v) < MIN_VEC_MAG) break;

        glm::vec2 np = pathStep(p, ft, h);
        float nt = ft + h;

        // stopping: boundary in space or time
        if (!insideGrid(np.x, np.y)) break;
        if (nt < 0.0f || nt > (float)(num_timesteps - 1)) break;

        float seg_len = glm::length(np - p);
        length += seg_len;
        if (length > MAX_LENGTH) break;

        glm::vec2 cur_ndc = gridToNDC(np.x, np.y);
        appendSegment(out, prev_ndc, cur_ndc);
        prev_ndc = cur_ndc;
        p = np;
        ft = nt;
    }
}

// ---------------------------------------------------------------------------
// Public entry points used by mouse / timestep callbacks
// ---------------------------------------------------------------------------

void recomputeAllStreamlines(void)
{
    streamline_vertices.clear();
    for (size_t i = 0; i < streamline_seeds.size(); ++i) {
        glm::vec2 p0(streamline_seeds[i].x, streamline_seeds[i].y);
        integrateStream(p0, loaded_timestep, +1.0f, streamline_vertices);
        integrateStream(p0, loaded_timestep, -1.0f, streamline_vertices);
    }
    uploadOverlay(streamVAO, streamVBO, streamline_vertices);
}

static void appendOnePathline(glm::vec2 p0, int t0)
{
    integratePath(p0, (float)t0, +1.0f, pathline_vertices);
    integratePath(p0, (float)t0, -1.0f, pathline_vertices);
}

void computeStreamline(int x, int y)
{
    if (!grid_data_loaded) return;
    glm::vec2 p0((float)x, (float)y);

    // store seed and append segments to the existing buffer
    Seed s; s.x = p0.x; s.y = p0.y; s.t = loaded_timestep;
    streamline_seeds.push_back(s);

    integrateStream(p0, loaded_timestep, +1.0f, streamline_vertices);
    integrateStream(p0, loaded_timestep, -1.0f, streamline_vertices);

    uploadOverlay(streamVAO, streamVBO, streamline_vertices);
}

void computePathline(int x, int y, int t)
{
    if (!grid_data_loaded) return;
    glm::vec2 p0((float)x, (float)y);

    Seed s; s.x = p0.x; s.y = p0.y; s.t = t;
    pathline_seeds.push_back(s);

    appendOnePathline(p0, t);
    uploadOverlay(pathVAO, pathVBO, pathline_vertices);
}

// ---------------------------------------------------------------------------
// Glyphs
// ---------------------------------------------------------------------------

// build one arrow (shaft + 2 head segments) in NDC into out
static void buildArrow(glm::vec2 grid_p, glm::vec2 v, float base_len,
                       std::vector<float>& out)
{
    float speed = glm::length(v);
    if (speed < MIN_VEC_MAG) return;

    glm::vec2 dir = v / speed;

    // length in grid units
    float len;
    if (en_constant_length) {
        len = base_len;
    } else {
        // crude scaling: clamp speed*0.3 within [0.2..1.0] of base_len
        float s = std::min(speed * 0.3f, 1.0f);
        s = std::max(s, 0.2f);
        len = base_len * s;
    }

    glm::vec2 tip   = grid_p + len * dir;
    glm::vec2 perp(-dir.y, dir.x);
    float head_back = len * 0.30f;
    float head_side = len * 0.18f;
    glm::vec2 h1 = tip - head_back * dir + head_side * perp;
    glm::vec2 h2 = tip - head_back * dir - head_side * perp;

    glm::vec2 a = gridToNDC(grid_p.x, grid_p.y);
    glm::vec2 b = gridToNDC(tip.x,    tip.y);
    glm::vec2 c = gridToNDC(h1.x,     h1.y);
    glm::vec2 d = gridToNDC(h2.x,     h2.y);
    appendSegment(out, a, b);
    appendSegment(out, b, c);
    appendSegment(out, b, d);
}

void drawGlyphs()
{
    if (!grid_data_loaded) return;

    glyph_vertices.clear();
    // sampling_rate is the stride between sampled cells
    int stride = std::max(1, sampling_rate);
    float base_len = (float)stride * 0.7f;

    for (int y = stride / 2; y < vol_dim[1]; y += stride) {
        for (int x = stride / 2; x < vol_dim[0]; x += stride) {
            glm::vec2 v = vectorAtCell(x, y, loaded_timestep);
            buildArrow(glm::vec2((float)x, (float)y), v, base_len, glyph_vertices);
        }
    }

    uploadOverlay(glyphVAO, glyphVBO, glyph_vertices);
    drawOverlay(glyphVAO, (GLsizei)(glyph_vertices.size() / 2),
                glm::vec4(1.0f, 1.0f, 1.0f, 1.0f));
}

// ---------------------------------------------------------------------------
// Misc helpers
// ---------------------------------------------------------------------------

static void clearAllSeeds(void)
{
    streamline_seeds.clear();
    pathline_seeds.clear();
    streamline_vertices.clear();
    pathline_vertices.clear();
    uploadOverlay(streamVAO, streamVBO, streamline_vertices);
    uploadOverlay(pathVAO,   pathVBO,   pathline_vertices);
}

static void seedRakeAt(float gx, float gy)
{
    // emit N seeds along a line perpendicular to the click axis spanning the
    // full grid extent. Vertical rake: fixed x = gx, varying y. Horizontal: opposite.
    const int N = 12;
    if (rake_mode == RAKE_VERTICAL) {
        float spacing = (float)(vol_dim[1] - 1) / (float)(N + 1);
        for (int i = 1; i <= N; ++i) {
            float yy = spacing * (float)i;
            if (en_streamline) computeStreamline((int)gx, (int)yy);
            if (en_pathline)   computePathline((int)gx, (int)yy, loaded_timestep);
        }
    } else if (rake_mode == RAKE_HORIZONTAL) {
        float spacing = (float)(vol_dim[0] - 1) / (float)(N + 1);
        for (int i = 1; i <= N; ++i) {
            float xx = spacing * (float)i;
            if (en_streamline) computeStreamline((int)xx, (int)gy);
            if (en_pathline)   computePathline((int)xx, (int)gy, loaded_timestep);
        }
    }
}

// cycle clear colors
static void nextClearColor()
{
    clearColor = (++clearColor) % 3;

    switch(clearColor)
    {
        case 0:
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            break;
        case 1:
            glClearColor(0.2f, 0.2f, 0.3f, 1.0f);
            break;
        default:
            glClearColor(0.7f, 0.7f, 0.7f, 1.0f);
            break;
    }
}


// callbacks
// framebuffer to fix viewport
void frameBufferCallback(GLFWwindow* window, int width, int height)
{
    view_width = width;
    view_height = height;
    glViewport(0, 0, width, height);
}

// key callback to take user inputs for both windows
void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
    if (action != GLFW_RELEASE) {
        const char* status[ 2 ];
        status[ 0 ] = "disabled";
        status[ 1 ] = "enabled";

        switch (key) {
            case '1':
                toggle_xy = 0;
                LoadData( filenames[ 0 ] );
                loaded_file = 0;
                fprintf( stderr, "Loading " );
                fprintf( stderr, "%s", filenames[ 0 ] );
                fprintf( stderr, " dataset.\n");
                clearAllSeeds();
                break;
            case '2':
                toggle_xy = 0;
                LoadData(filenames[ 1 ] );
                loaded_file = 1;
                fprintf( stderr, "Loading " );
                fprintf( stderr, "%s", filenames[ 1 ] );
                fprintf( stderr, " dataset.\n");
                clearAllSeeds();
                break;
            case '3':
                toggle_xy = 1;
                LoadData( filenames[ 2 ] );
                loaded_file = 2;
                fprintf( stderr, "Loading " );
                fprintf( stderr, "%s", filenames[ 2 ] );
                fprintf( stderr, " dataset.\n");
                clearAllSeeds();
                break;
            case '0':
                if( num_timesteps > 1 ){
                    loadNextTimestep();
                    fprintf( stderr, "Timestep %d.\n", loaded_timestep );
                    // streamlines depend on the current slice -> recompute
                    recomputeAllStreamlines();
                }
                break;
            case GLFW_KEY_A:
                en_arrow = !en_arrow;
                fprintf(stderr, "%s drawing arrows.\n", en_arrow? "enabling" : "disabling");
                break;
            case GLFW_KEY_S:
                current_scalar_field = (current_scalar_field + 1)%num_scalar_fields;
                DownloadScalarFieldAsTexture();
                fprintf( stderr, "Scalar field changed.\n");
                break;
            case GLFW_KEY_B:
                nextClearColor();
                fprintf( stderr, "Next clear color.\n");
                break;
            case GLFW_KEY_EQUAL:
                sampling_rate = std::min(sampling_rate + 5, 100);
                fprintf(stderr, "Increasing sampling rate to %d.\n", sampling_rate);
                break;
            case GLFW_KEY_MINUS:
                sampling_rate = std::max(sampling_rate - 5, 5);
                fprintf(stderr, "Decreasing sampling rate to: %d.\n", sampling_rate);
                break;
            case GLFW_KEY_I:
                dt = std::min(dt + 0.005f, 1.0f);
                fprintf(stderr, "Increase dt: %.3f\n", dt);
                break;
            case GLFW_KEY_K:
                dt = std::max(dt - 0.005f, 0.0001f);
                fprintf(stderr, "Decrease dt: %.3f\n", dt);
                break;
            case GLFW_KEY_T:
                en_streamline = !en_streamline;
                fprintf(stderr, "%s drawing streamlines.\n", en_streamline? "enabling" : "disabling");
                break;
            case GLFW_KEY_P:
                en_pathline = !en_pathline;
                fprintf(stderr, "%s drawing pathlines.\n", en_pathline? "enabling" : "disabling");
                break;
            case GLFW_KEY_C:
                colormap_mode = (colormap_mode + 1) % 3;
                fprintf(stderr, "Colormap: %s\n",
                        colormap_mode == CMAP_OFF ? "off (grayscale)" :
                        colormap_mode == CMAP_RAINBOW ? "rainbow" : "cool-warm");
                break;
            case GLFW_KEY_LEFT_BRACKET:
                blend_factor = std::max(0.0f, blend_factor - 0.1f);
                fprintf(stderr, "Blend factor: %.2f\n", blend_factor);
                break;
            case GLFW_KEY_RIGHT_BRACKET:
                blend_factor = std::min(1.0f, blend_factor + 0.1f);
                fprintf(stderr, "Blend factor: %.2f\n", blend_factor);
                break;
            case GLFW_KEY_M:
                integration_method = (integration_method + 1) % 3;
                fprintf(stderr, "Integration: %s\n",
                        integration_method == INTEG_EULER ? "Euler" :
                        integration_method == INTEG_RK2   ? "RK2"   : "RK4");
                // streamlines depend on integration method; recompute
                recomputeAllStreamlines();
                break;
            case GLFW_KEY_L:
                en_constant_length = !en_constant_length;
                fprintf(stderr, "Arrow length: %s\n",
                        en_constant_length ? "constant" : "speed-scaled");
                break;
            case GLFW_KEY_R:
                rake_mode = (rake_mode + 1) % 3;
                fprintf(stderr, "Rake mode: %s\n",
                        rake_mode == RAKE_OFF        ? "off" :
                        rake_mode == RAKE_VERTICAL   ? "vertical"
                                                     : "horizontal");
                break;
            case GLFW_KEY_X:
                clearAllSeeds();
                fprintf(stderr, "Cleared all seeds.\n");
                break;
            case GLFW_KEY_Q:
            case GLFW_KEY_ESCAPE:
                exit( 0 );
                break;
            default:
                fprintf( stderr, "\nKeyboard commands:\n\n"
                                 "1, load %s dataset\n"
                                 "2, load %s dataset\n"
                                 "3, load %s dataset\n"
                                 "0, cycle through timesteps\n"
                                 "b, switch background color\n"
                                 "a, en-/disable arrows.\n"
                                 "l, toggle arrow length (constant/speed)\n"
                                 "t, en-/disable streamlines.\n"
                                 "p, en-/disable pathlines.\n"
                                 "m, cycle integration method (Euler/RK2/RK4)\n"
                                 "r, cycle rake mode (off/vertical/horizontal)\n"
                                 "x, clear all seeds\n"
                                 "s, switch scalar field\n"
                                 "c, cycle colormap (off/rainbow/cool-warm)\n"
                                 "[ ], decrease/increase blend factor\n"
                                 "+, increase sampling rate.\n"
                                 "-, decrease sampling rate.\n"
                                 "i, increase dt.\n"
                                 "k, decrease dt.\n"
                                 "left-click, seed streamline/pathline\n"
                                 "q, <esc> - Quit\n",
                         filenames[0], filenames[1], filenames[2] );
                break;
        }
    }
}

// mouse callback to seed streamlines/pathlines
void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
    {
        if (!grid_data_loaded) return;

        double xpos, ypos;
        glfwGetCursorPos(window, &xpos, &ypos);

        int ww, hh;
        glfwGetWindowSize(window, &ww, &hh);
        if (ww <= 0 || hh <= 0) return;

        // convert screen coords -> grid coords (y axis flipped)
        float gx = (float)xpos / (float)ww * (float)(vol_dim[0] - 1);
        float gy = (1.0f - (float)ypos / (float)hh) * (float)(vol_dim[1] - 1);

        if (!insideGrid(gx, gy)) return;

        if (rake_mode != RAKE_OFF) {
            seedRakeAt(gx, gy);
        } else {
            if (en_streamline) computeStreamline((int)gx, (int)gy);
            if (en_pathline)   computePathline((int)gx, (int)gy, loaded_timestep);
        }
    }
}

// glfw error callback
static void errorCallback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

// data

void loadNextTimestep( void )
{
    loaded_timestep = ( loaded_timestep + 1 ) % num_timesteps;
    DownloadScalarFieldAsTexture();
}


/*
 * load .gri dataset
 * This only reads the header information and calls the dat loader
 * For now we ignore the grid data and assume a rectangular grid
 */
void LoadData( char* base_filename )
{
    //reset
    reset_rendering_props();

    char filename[ 80 ];
    strcpy( filename, base_filename );
    strcat( filename, ".gri");

    fprintf( stderr, "loading grid file %s\n", filename );

    // open grid file, read only, binary mode
    FILE* fp = fopen( filename, "rb" );
    if ( fp == NULL ) {
        fprintf( stderr, "Cannot open file %s for reading.\n", filename );
        return;
    }

    // read header
    char header[ 40 ];
    fread( header, sizeof( char ), 40, fp );
    sscanf( header, "SN4DB %d %d %d %d %d %f",
            &vol_dim[ 0 ], &vol_dim[ 1 ], &vol_dim[ 2 ],
            &num_scalar_fields, &num_timesteps ,&timestep );

    fprintf( stderr, "dimensions: x: %d, y: %d, z: %d.\n", vol_dim[ 0 ], vol_dim[ 1 ], vol_dim[ 2 ] );
    fprintf( stderr, "additional info: # scalar fields: %d, # timesteps: %d, timestep: %f.\n", num_scalar_fields, num_timesteps, timestep );

    // read data
    char dat_filename[ 80 ];
    strcpy( dat_filename, base_filename );

    if( num_timesteps <= 1 ){

        strcat( dat_filename, ".dat");

    } else {

        strcat( dat_filename, ".00000.dat");

    }

    loaded_timestep = 0;
    LoadVectorData( base_filename );

    glfwSetWindowSize(window, vol_dim[ 0 ], vol_dim[ 1 ] );
    grid_data_loaded = true;
}

/*
 * load .dat dataset
 * loads vector and scalar fields
 */
void LoadVectorData( const char* filename )
{
    fprintf( stderr, "loading scalar file %s\n", filename );

    // open data file, read only, binary mode
    char ts_name[ 80 ];
    if( num_timesteps > 1 )
    {
        sprintf( ts_name, "%s.%.5d.dat", filename, 0 );
    }
    else
        sprintf( ts_name, "%s.dat",filename);

    FILE* fp = fopen( ts_name, "rb" );
    if ( fp == NULL ) {
        fprintf( stderr, "Cannot open file %s for reading.\n", filename );
        return;
    }
    else
    {
        fclose( fp );
    }

    data_size = vol_dim[ 0 ] * vol_dim[ 1 ] * vol_dim[ 2 ];

    if (!vector_array) {
        delete[] vector_array;
        vector_array = NULL;
    }
    // dim.xyz * vector.xyz * timesteps
    vector_array = new float[ data_size * 3 * num_timesteps];

    // read data
    if (!scalar_fields) {
        delete[] scalar_fields;
        scalar_fields = NULL;
        delete[] scalar_bounds;
        scalar_bounds = NULL;
    }
    // dim.xyz * scalarfields(2) * timesteps
    scalar_fields = new float[ data_size * num_scalar_fields * num_timesteps ];
    scalar_bounds = new float[ 2 * num_scalar_fields * num_timesteps ];

    int num_total_fields = num_scalar_fields + 3; // scalar fields + vec.xyz
    float *tmp = new float[ data_size * num_total_fields * num_timesteps ];

    for( int k = 0 ; k < num_timesteps; k++ )
    {
        char times_name[ 80 ];
        if( num_timesteps > 1 )
        {
            sprintf( times_name, "%s.%.5d.dat", filename, k );
            fp = fopen( times_name, "rb" );
        }
        else
            fp = fopen( ts_name, "rb" );
        // read scalar data
        fread( &tmp[k*data_size*num_total_fields], sizeof( float ), ( data_size * num_total_fields ), fp );

        // close file
        fclose( fp );

        // copy and scan for min and max values
        for( int  i = 0; i < num_scalar_fields; i++ ){

            float min_val = 99999.9f;
            float max_val = 0.0f;

            float avg = 0.0;

            int offset = i * data_size * num_timesteps;

            for( int j = 0; j < data_size; j++ ){

                float val = tmp[ j * num_total_fields + 3 + i + k*data_size*num_total_fields ];

                scalar_fields[ j + k*data_size + offset ] = val;

                if( toggle_xy ) {
                    // overwrite
                    if( i == 0 ){
                        vector_array[ 3*j + 0 + 3*k*data_size ] = tmp[ j * num_total_fields + 1 + k*data_size*num_total_fields ];//toggle x and y components in vector field
                        vector_array[ 3*j + 1 + 3*k*data_size ] = tmp[ j * num_total_fields + 0 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 2 + 3*k*data_size ] = tmp[ j * num_total_fields + 2 + k*data_size*num_total_fields ];
                    }
                } else {
                    // overwrite
                    if( i == 0 ){
                        vector_array[ 3*j + 0 + 3*k*data_size ] = tmp[ j * num_total_fields + 0 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 1 + 3*k*data_size ] = tmp[ j * num_total_fields + 1 + k*data_size*num_total_fields ];
                        vector_array[ 3*j + 2 + 3*k*data_size ] = tmp[ j * num_total_fields + 2 + k*data_size*num_total_fields ];
                    }
                }

                min_val = std::min( val, min_val );
                max_val = std::max( val, max_val );

                avg += scalar_fields[ offset + j + k*data_size ] / data_size;
            }
            scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = min_val;
            scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = max_val;
        }

        // normalize
        for( int  i = 0; i < num_scalar_fields; i++ ){

            int offset = i * data_size * num_timesteps;

            float lower_bound = scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ];
            float upper_bound = scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ];

            // scale between [0..1] where 1 is original zero
            // the boundary of the bigger abs border will be used to scale
            // meaning one boundary will likely not be hit i.e real scale
            // for -50..100 will be [0.25..1.0]
            if( lower_bound < 0.0 && upper_bound > 0.0 ){

                float scale = 0.5f / std::max( -lower_bound, upper_bound );

                for( int j = 0; j < data_size; j++ ){

                    scalar_fields[ offset + j + k*data_size ] = 0.5f + scalar_fields[ offset + j + k*data_size ] * scale;
                }
                scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = 0.5f + scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] * scale;
                scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = 0.5f + scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] * scale;


                // scale between [0..1]
            } else {

                float sign = upper_bound <= 0.0 ? -1.0f : 1.0f;

                float scale = 1.0f / ( upper_bound - lower_bound ) * sign;

                for( int j = 0; j < data_size; j++ ){

                    scalar_fields[ offset + j + k*data_size ] = ( scalar_fields[ offset + j + k*data_size ] - lower_bound ) * scale;
                }
                scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] = ( scalar_bounds[ 2 * i     + k*num_scalar_fields*2 ] + lower_bound ) * scale; //should be 0.0
                scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] = ( scalar_bounds[ 2 * i + 1 + k*num_scalar_fields*2 ] + lower_bound ) * scale; //should be 1.0
            }
        }
    }
    delete[] tmp;
    DownloadScalarFieldAsTexture();

    scalar_data_loaded = true;
}


void DownloadScalarFieldAsTexture( void )
{
    fprintf( stderr, "downloading scalar field to 2D texture\n" );

    // (re)generate the 2D texture each time we change slice/scalar field
    glGenTextures( 1, &scalar_field_texture );
    glBindTexture( GL_TEXTURE_2D, scalar_field_texture );

    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
    glTexParameteri( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );

    int datasize = vol_dim[0] * vol_dim[1];
    // single-channel float texture; the fragment shader reads .r
    glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
    glTexImage2D( GL_TEXTURE_2D, 0, GL_R16, vol_dim[ 0 ], vol_dim[ 1 ], 0,
                  GL_RED, GL_FLOAT,
                  &scalar_fields[ (loaded_timestep + current_scalar_field * num_timesteps)*datasize ] );
}

char *textFileRead( char *fn ){

    FILE *fp;
    char *content = NULL;

    int count=0;

    if (!fn) {
        fp = fopen(fn,"rt");

        if (!fp) {

            fseek(fp, 0, SEEK_END);
            count = ftell(fp);
            rewind(fp);

            if (count > 0) {
                content = (char *)malloc(sizeof(char) * (count+1));
                count = fread(content,sizeof(char),count,fp);
                content[count] = '\0';
            }
            fclose(fp);
        }
    }
    return content;
}


// initializations
// init application
bool initApplication(int argc, char **argv)
{

    std::string version((const char *)glGetString(GL_VERSION));
    std::stringstream stream(version);
    unsigned major, minor;
    char dot;

    stream >> major >> dot >> minor;

    assert(dot == '.');
    if (major > 3 || (major == 2 && minor >= 0)) {
        std::cout << "OpenGL Version " << major << "." << minor << std::endl;
    } else {
        std::cout << "The minimum required OpenGL version is not supported on this machine. Supported is only " << major << "." << minor << std::endl;
        return false;
    }

    return true;
}

void reset_rendering_props( void )
{
    num_scalar_fields = 0;
}

// set up the scene
void setup() {
    LoadData( filenames[ 0 ] );
    loaded_file = 0;

    DownloadScalarFieldAsTexture();


    // compile & link shader
    std::string vs = resolveAsset("shaders/vertex.vs");
    std::string fs = resolveAsset("shaders/fragment.fs");
    vectorProgram.compileShader(vs.c_str());
    vectorProgram.compileShader(fs.c_str());
    vectorProgram.link();

    // make quad to render texture
    // see: vboquad.h and vboquad.cpp
    quad.init();

    // overlay buffers (glyphs, streamlines, pathlines)
    initOverlay(glyphVAO,  glyphVBO);
    initOverlay(streamVAO, streamVBO);
    initOverlay(pathVAO,   pathVBO);
}

// rendering
void render() {
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

    // draw the texture (background scalar field, with optional colormap)
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, scalar_field_texture);
    vectorProgram.use();

    model = glm::mat4(1);

    vectorProgram.setUniform("vertexColor", glm::vec4(0));
    vectorProgram.setUniform("model", model);
    vectorProgram.setUniform("colormapMode", colormap_mode);
    vectorProgram.setUniform("blendFactor",  blend_factor);
    vectorProgram.setUniform("txtr",         0);

    quad.render();

    // overlays use a solid color, no texture
    if (en_arrow)      drawGlyphs();
    if (en_streamline) drawOverlay(streamVAO,
                                   (GLsizei)(streamline_vertices.size() / 2),
                                   glm::vec4(1.0f, 1.0f, 0.0f, 1.0f));
    if (en_pathline)   drawOverlay(pathVAO,
                                   (GLsizei)(pathline_vertices.size() / 2),
                                   glm::vec4(1.0f, 0.2f, 0.2f, 1.0f));
}

// entry point
int main(int argc, char** argv)
{
    // init variables
    view_width = 0;
    view_height = 0;

    toggle_xy = 0;

    en_arrow = false;
    en_streamline = false;
    en_pathline = false;
    sampling_rate = 15;
    dt = 0.1f;

    integration_method = INTEG_EULER;
    en_constant_length = true;
    rake_mode          = RAKE_OFF;

    colormap_mode = CMAP_OFF;
    blend_factor  = 1.0f;

    reset_rendering_props();

    vector_array = NULL;
    scalar_fields = NULL;
    scalar_bounds = NULL;
    grid_data_loaded = false;
    scalar_data_loaded = false;
    current_scalar_field = 0;
    clearColor = 0;


    // resolve dataset paths relative to wherever the executable was launched
    static std::string f0 = resolveAsset("data/block/c_block");
    static std::string f1 = resolveAsset("data/tube/tube");
    static std::string f2 = resolveAsset("data/hurricane/hurricane_p_tc");
    filenames[ 0 ] = (char*)f0.c_str();
    filenames[ 1 ] = (char*)f1.c_str();
    filenames[ 2 ] = (char*)f2.c_str();



    // set glfw error callback
    glfwSetErrorCallback(errorCallback);

    // init glfw
    if (!glfwInit()) {
        exit(EXIT_FAILURE);
    }

    // request an OpenGL 4.1 core context (max supported on macOS; works elsewhere too)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE,        GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    // init glfw window
    window = glfwCreateWindow(gWindowWidth, gWindowHeight, "AMCS/CS247 Scientific Visualization", nullptr, nullptr);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // set GLFW callback functions
    glfwSetKeyCallback(window, keyCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);
    glfwSetFramebufferSizeCallback(window, frameBufferCallback);

    // make context current (once is sufficient)
    glfwMakeContextCurrent(window);

    // get the frame buffer size
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    // init the OpenGL API (we need to do this once before any calls to the OpenGL API)
    gladLoadGL();

    // init our application
    if (!initApplication(argc, argv)) {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }


    // set up the scene
    setup();

    // print menu
    keyCallback(window, GLFW_KEY_BACKSLASH, 0, GLFW_PRESS, 0);

    // start traversing the main loop
    // loop until the user closes the window
    while (!glfwWindowShouldClose(window))
    {
        // clear frame buffer
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // render one frame
        render();

        // swap front and back buffers
        glfwSwapBuffers(window);

        // poll and process input events (keyboard, mouse, window, ...)
        glfwPollEvents();
    }

    glfwTerminate();
    return EXIT_SUCCESS;
}
