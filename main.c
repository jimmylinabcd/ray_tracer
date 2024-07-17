#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>

#include <cglm/cglm.h>
#include <cglm/call.h>

#define GUILITE_ON
#include "GuiLite.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

typedef vec3                    mat3[3];

#ifdef __AVX__
typedef CGLM_ALIGN_IF(32) vec4  mat4[4];
#else
typedef CGLM_ALIGN_IF(16) vec4  mat4[4];
#endif


#define OUTPUT_IMAGE_HEIGHT 1080
#define OUTPUT_IMAGE_WIDTH 1920

#define AA_SAMPLES 1
#define BOUNCE_DEPTH 60

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t);


void button_callback(void *ctx, int id)
{
    printf("Button clicked ID: %d\n", id);
}

int main(){
    // Initialize the GUI
    gui_t gui;
    gui_init(&gui, 800, 600, "Test");

    // Create a button
    gui_element_t button = gui_create_button(&gui, 350, 250, 100, 50, "click", button_callback, NULL);

    // Main loop
    while (gui_process(&gui))
    {
        // Clear the screen
        gui_clear(&gui);

        // Render the button
        gui_render_element(&gui, &button);

        // Update the display
        gui_display(&gui);
    }

    // Clean up
    gui_destroy(&gui);
    return 0;

    printf("Running Ray Tracer\n");

    srand(time(NULL));

    // Camera
    double focal_length = 1.0;

    double viewport_height = 2.0;
    double viewport_width = (OUTPUT_IMAGE_WIDTH/OUTPUT_IMAGE_HEIGHT) * viewport_height;

    vec3 origin = {0, 0, 0};
    vec3 horizontal =  {viewport_width, 0, 0};
    vec3 verticle = {0, viewport_height, 0};

    vec3 lower_left_corner;
    lower_left_corner[0] = origin[0] - horizontal[0]/2 - verticle[0]/2;
    lower_left_corner[1] = origin[1] - horizontal[1]/2 - verticle[1]/2;
    lower_left_corner[2] = origin[2] - horizontal[2]/2 - verticle[2]/2 - focal_length;
    // Scene Description

    // Loading Object

    unsigned char buffer[4];

    char file_name[] = "knight.stl";

    FILE *ptr;

    ptr = fopen(file_name,"rb");

    // Jump Past 80 Byte Header
    fseek(ptr, 80, SEEK_SET);

    // Read Number of Triangles
    unsigned int number_triangles;
    fread(&number_triangles,sizeof(unsigned int),1,ptr);
  
    printf("Loading %d triangles from %s...\n", number_triangles, file_name);

    vec3 *vertices = (vec3  *)malloc(number_triangles  * 3 * sizeof(vec3));

    // Load vertecies into memeory

    // bounding box
    float maxX = 0;
	float maxY = 0;
	float minX = 0;
	float minY = 0;

    for (unsigned int i = 0; i < number_triangles; i++) {
        // discard normal vector 12 bytes
        // Maybe we can use this?
        fseek(ptr, 12, SEEK_CUR);

        // Read vertices 3 * 12 bytes
        fread(&vertices[i * 3], sizeof(vec3), 3, ptr);

		if(vertices[i * 3][0] > maxX){
			maxX = vertices[i * 3][0];
		}else if(vertices[i * 3][0] < minX){
			minX = vertices[i * 3][0];
		}


		if(vertices[i * 3][1] > maxY){
			maxY = vertices[i * 3][1];
		}else if(vertices[i * 3][1] < minY){
			minX = vertices[i * 3][1];
		}

        // Discard attribute 2 byte 
        fseek(ptr, 2, SEEK_CUR);
    }

    fclose(ptr);

    // Build save or maybe load bvh?

    // Image

    char* frameData = malloc(OUTPUT_IMAGE_WIDTH * OUTPUT_IMAGE_HEIGHT * 3 * sizeof(char));

    unsigned int index = 0;

    for (unsigned int y = OUTPUT_IMAGE_HEIGHT - 1; y >=0 ; y--) {
        for (unsigned int x = 0; x < OUTPUT_IMAGE_WIDTH; x++) {

            vec3 pixel_colour = {0, 0, 0};

            for (int i = 0; i < AA_SAMPLES; i++){

                double u = (double) x + ((rand()/ (RAND_MAX + 1.0))) / (OUTPUT_IMAGE_WIDTH - 1);
                double v = (double) y + ((rand()/ (RAND_MAX + 1.0))) / (OUTPUT_IMAGE_HEIGHT - 1);
                

                vec3 ray_direction; 
                ray_direction[0] = lower_left_corner[0] + u * horizontal[0] + v * verticle[0] - origin[0];
                ray_direction[1] = lower_left_corner[1] + u * horizontal[1] + v * verticle[1] - origin[1];
                ray_direction[2] = lower_left_corner[2] + u * horizontal[2] + v * verticle[2] - origin[2];
                
                vec3 vertex1 = {-0.5, -0.5, -1};
                vec3 vertex2 = {0.5, -0.5, -1};
                vec3 vertex3 = {0, 0.5, -1};

                float distance;

                glm_vec3_norm2(ray_direction);

                ray_triangle_intersection(origin, ray_direction, vertex1, vertex2, vertex3, &distance);

                if (distance){
                    pixel_colour[0] = 255;
                }
            }

            frameData[index] = pixel_colour[0];
            frameData[index + 1] = pixel_colour[1];
            frameData[index + 2] = pixel_colour[2];

            index += 3*sizeof(char);
        }
    }

    stbi_write_png("output.png", OUTPUT_IMAGE_WIDTH, OUTPUT_IMAGE_HEIGHT, 3, frameData, 0);

    return 0;
}

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t) {
    const float epsilon = 0.000001;

    vec3 edge1, edge2, ray_cross_e2, s, s_cross_e1;

    glm_vec3_sub(vertex2, vertex1, edge1);
    glm_vec3_sub(vertex3, vertex1, edge2);

    glm_vec3_cross(ray_vector, edge2, ray_cross_e2);

    float det = glm_vec3_dot(edge1, ray_cross_e2);

    if (det > -epsilon && det < epsilon) {
        // Ray is parallel to triangle
        *t = -1.0f; // Indicate no intersection
        return;
    }

    float inv_det = 1.0 / det;
    glm_vec3_sub(ray_origin, vertex1, s);

    float u = inv_det * glm_vec3_dot(s, ray_cross_e2);

    if (u < 0 || u > 1) {
        *t = -1.0f; // Indicate no intersection
        return;
    }

    glm_vec3_cross(edge2, s, s_cross_e1);
    float v = inv_det * glm_vec3_dot(ray_vector, s_cross_e1);

    if (v < 0 || u + v > 1) {
        *t = -1.0f; // Indicate no intersection
        return;
    }

    *t = inv_det * glm_vec3_dot(edge1, s_cross_e1);

    if (*t > epsilon) {
        return;
    }
}