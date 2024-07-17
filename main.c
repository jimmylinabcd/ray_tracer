#include <stdio.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <cglm/cglm.h>
#include <cglm/call.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

typedef vec3 mat3[3];

#ifdef __AVX__
typedef CGLM_ALIGN_IF(32) vec4 mat4[4];
#else
typedef CGLM_ALIGN_IF(16) vec4 mat4[4];
#endif

#define OUTPUT_IMAGE_HEIGHT 1080
#define OUTPUT_IMAGE_WIDTH 1920

#define AA_SAMPLES 1
#define BOUNCE_DEPTH 60

typedef struct {
    vec3 min;
    vec3 max;
} bounding_box;

typedef struct BVH_node{
    bounding_box bdbox;
    struct BVH_node* left;
    struct BVH_node* right;
    unsigned int obj_start;
    unsigned int obj_end;
} BVH_node;

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t, vec3* out_intersection);
void ray_plane_intersection(vec3 ray_origin, vec3 ray_vector, vec3 plane_point, vec3 plane_normal, float* t, vec3* out_intersection);
void ray_box_intersection(vec3 ray_origin, vec3 ray_direction, vec3 box_min, vec3 box_max, float* tmin, float* tmax);
vec3 trace(vec3 ray_origin, vec3 ray_direction);

int main() {
    printf("Running Ray Tracer\n");

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


    srand(time(NULL));

    // Camera
    double focal_length = 1.0;

    double viewport_height = 2.0;
    double viewport_width = (OUTPUT_IMAGE_WIDTH / OUTPUT_IMAGE_HEIGHT) * viewport_height;

    vec3 origin = {0, 0, 0};
    vec3 horizontal = {viewport_width, 0, 0};
    vec3 vertical = {0, viewport_height, 0};

    vec3 lower_left_corner;
    glm_vec3_copy(origin, lower_left_corner);
    glm_vec3_sub(lower_left_corner, horizontal, lower_left_corner);
    glm_vec3_sub(lower_left_corner, vertical, lower_left_corner);
    lower_left_corner[2] -= focal_length;

    // Scene Description

    // Image
    unsigned char* frameData = malloc(OUTPUT_IMAGE_WIDTH * OUTPUT_IMAGE_HEIGHT * 3 * sizeof(char));

    unsigned int index = 0;

    for (int y = OUTPUT_IMAGE_HEIGHT - 1; y >= 0; y--) {
        for (unsigned int x = 0; x < OUTPUT_IMAGE_WIDTH; x++) {
            vec3 pixel_colour = {0, 0, 0};

            for (int i = 0; i < AA_SAMPLES; i++) {
                double u = (double) x / (OUTPUT_IMAGE_WIDTH - 1);
                double v = (double) y / (OUTPUT_IMAGE_HEIGHT - 1);

                vec3 ray_direction;
                glm_vec3_copy(lower_left_corner, ray_direction);
                vec3 scaled_horizontal, scaled_vertical;
                glm_vec3_scale(horizontal, u, scaled_horizontal);
                glm_vec3_scale(vertical, v, scaled_vertical);
                glm_vec3_add(ray_direction, scaled_horizontal, ray_direction);
                glm_vec3_add(ray_direction, scaled_vertical, ray_direction);
                glm_vec3_sub(ray_direction, origin, ray_direction);

                glm_vec3_normalize(ray_direction);

                vec3 temp = trace(origin, ray_direction);

                pixel_colour[0] = temp[0];
                pixel_colour[1] = temp[1];
                pixel_colour[2] = temp[2];
                
            }

            frameData[index] = (unsigned char) pixel_colour[0];
            frameData[index + 1] = (unsigned char) pixel_colour[1];
            frameData[index + 2] = (unsigned char) pixel_colour[2];

            index += 3;
        }
    }

    stbi_write_png("output.png", OUTPUT_IMAGE_WIDTH, OUTPUT_IMAGE_HEIGHT, 3, frameData, 0);

    free(frameData);

    return 0;
}

vec3 trace(vec3 ray_origin, vec3 ray_direction){
    vec3 vertex1 = {-0.5, -0.5, -1};
    vec3 vertex2 = {0.5, -0.5, -1};
    vec3 vertex3 = {0, 0.5, -1};

    float distance = -1.0f;
    vec3 intersection;
    ray_triangle_intersection(origin, ray_direction, vertex1, vertex2, vertex3, &distance, &intersection);
    if (distance > 0){
        vec3 pixel_colour = {255, 0, 0};
        return 
    }else{
        double a = 0.5*(ray[1] + 1.0);
        vec3 pixel_colour;
        // return (1.0-a)*color(1.0, 1.0, 1.0) + a*color(0.5, 0.7, 1.0);
        pixel_colour[0] = (1.0 - a) * 255 + a * 128;
        pixel_colour[1] = (1.0 - a) * 255 + a * 178;
        pixel_colour[2] = (1.0 - a) * 255 + a * 255;
        return pixel_colour;
    }
}

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t, vec3* out_intersection) {
    const float epsilon = 0.000001;

    vec3 edge1, edge2, ray_cross_e2, s, s_cross_e1;

    glm_vec3_sub(vertex2, vertex1, edge1);
    glm_vec3_sub(vertex3, vertex1, edge2);

    glm_vec3_cross(ray_vector, edge2, ray_cross_e2);

    float det = glm_vec3_dot(edge1, ray_cross_e2);

    if (det > -epsilon && det < epsilon) {
        *t = -1.0f;
        return;
    }

    float inv_det = 1.0 / det;
    glm_vec3_sub(ray_origin, vertex1, s);

    float u = inv_det * glm_vec3_dot(s, ray_cross_e2);

    if (u < 0 || u > 1) {
        *t = -1.0f;
        return;
    }

    glm_vec3_cross(s, edge1, s_cross_e1);
    float v = inv_det * glm_vec3_dot(ray_vector, s_cross_e1);

    if (v < 0 || u + v > 1) {
        *t = -1.0f;
        return;
    }

    *t = inv_det * glm_vec3_dot(edge2, s_cross_e1);

    if (*t > epsilon) {
        glm_vec3_scale(ray_vector, *t, *out_intersection);
        glm_vec3_add(ray_origin, *out_intersection, *out_intersection);
        return;
    }
}

void ray_plane_intersection(vec3 ray_origin, vec3 ray_vector, vec3 plane_point, vec3 plane_normal, float* t, vec3* out_intersection) {
    const float epsilon = 0.000001;
    float denom = glm_vec3_dot(plane_normal, ray_vector);

    // check if the ray is parallel to the plane
    if (fabs(denom) < epsilon) {
        *t = -1.0f;
        return;
    }

    vec3 p0l0;
    glm_vec3_sub(plane_point, ray_origin, p0l0);
    *t = glm_vec3_dot(p0l0, plane_normal) / denom;

    // check if the intersection is behind the ray origin
    if (*t < 0) {
        *t = -1.0f;
        return;
    }

    glm_vec3_scale(ray_vector, *t, *out_intersection);
    glm_vec3_add(ray_origin, *out_intersection, *out_intersection);
    return;
}

void ray_box_intersection(vec3 ray_origin, vec3 ray_direction, vec3 box_min, vec3 box_max, float* tmin, float* tmax) {
    vec3 inv_direction;
    glm_vec3_div(ray_direction, (vec3){1.0f, 1.0f, 1.0f}, inv_direction);

    vec3 tbot = glm_vec3_mul(glm_vec3_sub(box_min, ray_origin), inv_direction);
    vec3 ttop = glm_vec3_mul(glm_vec3_sub(box_max, ray_origin), inv_direction);

    vec3 tmin_v = glm_vec3_minv(ttop, tbot);
    vec3 tmax_v = glm_vec3_maxv(ttop, tbot);

    *tmin = glm_vec3_max3(tmin_v);
    *tmax = glm_vec3_min3(tmax_v);
}


