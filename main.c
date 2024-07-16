#include<stdio.h>
#include<stdbool.h>
#include <time.h>
#include <stdlib.h>

#include<cglm/cglm.h>
#include<cglm/call.h>

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

#define AA_SAMPLES 16
#define BOUNCE_DEPTH 60

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, vec3* target_ray);

int main(){
    printf("Running Ray Tracer\n");

    srand(time(NULL));

    // Camera
    double focal_length = 1.0;

    vec3 origin = {0, 0, 0};

    vec3 horizontal =  {OUTPUT_IMAGE_WIDTH, 0, 0};
    vec3 verticle = {0, OUTPUT_IMAGE_HEIGHT, 0};

    vec3 lower_left_corner = {origin[0] - horizontal[0]/2 - verticle[0]/2, origin[1] - horizontal[1]/2 - verticle[1]/2,  origin[2] - horizontal[2]/2 - verticle[2]/2 - focal_length};

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

    // Load vertecies into memeory

    fclose(ptr);

    // Build save or maybe load bvh?

    // Image

    char* frameData = malloc(OUTPUT_IMAGE_WIDTH * OUTPUT_IMAGE_HEIGHT * 3 * sizeof(char));

    unsigned int index = 0;

    for (int y = OUTPUT_IMAGE_HEIGHT - 1; y >=0 ; y--) {
        for (int x = 0; x < OUTPUT_IMAGE_WIDTH; x++) {

            vec3 pixel_colour;

            for (int i = 0; i < AA_SAMPLES; i++){

                double u = (double) x + ((rand()/ (RAND_MAX + 1.0))) / (OUTPUT_IMAGE_WIDTH - 1);
                double v = (double) y + ((rand()/ (RAND_MAX + 1.0))) / (OUTPUT_IMAGE_HEIGHT - 1);
                

                vec3 ray_direction; 
                ray_direction[0] = lower_left_corner[0] + u * horizontal[0] + v * verticle[0] - origin[0];
                ray_direction[1] = lower_left_corner[1] + u * horizontal[1] + v * verticle[1] - origin[1];
                ray_direction[2] = lower_left_corner[2] + u * horizontal[2] + v * verticle[2] - origin[2];
                
                vec3 vertex1 = {0, 0, -1};
                vec3 vertex2 = {1, 1, -1};
                vec3 vertex3 = {2, 2, -1};

                //float distance = 0;

                //if (glm_ray_triangle(origin, ray_direction, vertex1, vertex2, vertex3, &distance) == true){
                    //pixel_colour[0] = 255;
                //}

                vec3 unit_direction = unit_vector(r.direction());
                auto a = 0.5*(unit_direction.y() + 1.0);
                return (1.0-a)*color(1.0, 1.0, 1.0) + a*color(0.5, 0.7, 1.0);



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

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, vec3* target_ray) {
    const float epsilon = 0.000001;

    vec3 edge1, edge2, ray_cross_e2, s, s_cross_e1;

    glm_vec3_sub(vertex2, vertex1, edge1);
    glm_vec3_sub(vertex3, vertex1, edge2);

    glm_vec3_cross(ray_vector, edge2, ray_cross_e2);

    float det = glm_vec3_dot(edge1, ray_cross_e2);

    if (det > -epsilon && det < epsilon) {
        // Ray is parallel to triangle
        return;
    }

    float inv_det = 1.0 / det;
    glm_vec3_sub(ray_origin, vertex1, s);

    float u = inv_det * glm_vec3_dot(s, ray_cross_e2);

    if (u < 0 || u > 1) {
        return;
    }

    glm_vec3_cross(edge2, s, s_cross_e1);
    float v = inv_det * glm_vec3_dot(ray_vector, s_cross_e1);

    if (v < 0 || u + v > 1) {
        return;
    }

    float t = inv_det * glm_vec3_dot(edge1, s_cross_e1);

    if (t > epsilon) {
        vec3 temp;
        glm_vec3_scale(ray_vector, t, temp);
        glm_vec3_add(ray_origin, temp, temp);
        glm_vec3_copy(*target_ray, temp);
        return;
    }
}