#include <stdio.h>
#include <stdbool.h>
#include <string.h>
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

#define AA_SAMPLES 2
#define BOUNCE_DEPTH 1

typedef struct {
    vec3 min;
    vec3 max;
} bounding_box;

typedef struct BVH_node{
    bounding_box bdbox;
    struct BVH_node* left;
    struct BVH_node* right;
    // obj list needs to be sorted
    unsigned int triangle_index_start;
    unsigned int triangle_number;
} BVH_node;

typedef struct {
    vec3 vertex1;
    vec3 vertex2;
    vec3 vertex3;
    vec3 normal;
    vec3 colour;
} triangle;

typedef struct {
    unsigned int id;
    triangle* triangle_list;
    unsigned int list_size;
    BVH_node root;
} object;

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t, vec3* out_intersection);
void ray_plane_intersection(vec3 ray_origin, vec3 ray_vector, vec3 plane_point, vec3 plane_normal, float* t, vec3* out_intersection);
void ray_box_intersection(vec3 ray_origin, vec3 ray_direction, vec3 box_min, vec3 box_max, float* tmin, float* tmax, vec3* out_intersection);
void trace(vec3 ray_origin, vec3 ray_direction, vec3* out_pixel_colour, unsigned int depth, object* scene);
void compute_surface_normal(vec3 A, vec3 B, vec3 C, vec3 normal);
double randomRange(int min, int max);
double randomDouble();

void random_unit_vector(vec3 dest);

int main() {
    printf("Running Ray Tracer\n");

    // Scene Description

    object scene[2];

    vec3 origin = {0.0f, 0.0f, 0.0f};
    double focal_length = 1.0;
    double viewport_height = 2.0;
    double viewport_width = (OUTPUT_IMAGE_WIDTH / (double)OUTPUT_IMAGE_HEIGHT) * viewport_height;

    vec3 horizontal = {viewport_width, 0, 0};
    vec3 vertical = {0, viewport_height, 0};

    vec3 lower_left_corner;
    glm_vec3_copy(origin, lower_left_corner);
    glm_vec3_sub(lower_left_corner, horizontal, lower_left_corner);
    glm_vec3_sub(lower_left_corner, vertical, lower_left_corner);
    lower_left_corner[2] -= focal_length;

    glm_vec3_scale(lower_left_corner, 0.5f, lower_left_corner);
    glm_vec3_add(lower_left_corner, origin, lower_left_corner);


    // Building Scene

    // create test triangle
    triangle test_triangle;
    glm_vec3_copy((vec3){-0.5f, -0.5f, -1.0f}, test_triangle.vertex1);
    glm_vec3_copy((vec3){0.5f, -0.5f, -1.0f}, test_triangle.vertex2);
    glm_vec3_copy((vec3){0.0f, 0.5f, -1.0f}, test_triangle.vertex3);
    // calculate the normal
    vec3 temp_normal;
    compute_surface_normal(test_triangle.vertex1, test_triangle.vertex2, test_triangle.vertex3, temp_normal);


    // ground plane
    triangle test_triangle2;
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -1.0f}, test_triangle2.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -1.0f}, test_triangle2.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, test_triangle2.vertex3);
    vec3 temp_normal2;
    compute_surface_normal(test_triangle2.vertex1, test_triangle2.vertex2, test_triangle2.vertex3, temp_normal2);

    triangle test_triangle3;
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -1.0f}, test_triangle3.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, test_triangle3.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, test_triangle3.vertex3);
    vec3 temp_normal3;
    compute_surface_normal(test_triangle3.vertex1, test_triangle3.vertex2, test_triangle3.vertex3, temp_normal3);

    memcpy(test_triangle.normal, temp_normal, sizeof(vec3));
    memcpy(test_triangle2.normal, temp_normal2, sizeof(vec3));
    memcpy(test_triangle3.normal, temp_normal3, sizeof(vec3));
    
    glm_vec3_copy((vec3){255, 255, 255}, test_triangle.colour);
    glm_vec3_copy((vec3){255, 255, 255}, test_triangle2.colour);
    glm_vec3_copy((vec3){255, 255, 255}, test_triangle3.colour);

    // create test object
    scene[0].id = 0;
    //scene[0].triangle_list = &test_triangle;
    scene[0].triangle_list = malloc(sizeof(triangle) * 3);
    scene[0].triangle_list[0] = test_triangle;
    scene[0].triangle_list[1] = test_triangle2;
    scene[0].triangle_list[2] = test_triangle3;
    scene[0].list_size = 3;
    BVH_node test_bvh;
    
    test_bvh.left = NULL;
    test_bvh.right = NULL;
    test_bvh.triangle_index_start = 0;
    test_bvh.triangle_number = 1;
    
    bounding_box test_bounding_box;
    glm_vec3_copy((vec3){-0.5f, -0.5f, -1.0f}, test_bounding_box.min);
    glm_vec3_copy((vec3){0.5f, 0.5f, -1.0f}, test_bounding_box.max);
    test_bvh.bdbox = test_bounding_box;


    scene[0].root = test_bvh;


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

    scene[1].id = 1;
    scene[1].triangle_list = malloc(sizeof(triangle) * number_triangles);
    scene[1].list_size = number_triangles;
    scene[1].root = test_bvh;

    for (unsigned int i = 0; i < number_triangles; i++) {
        triangle temp;
        vec3 buffer[4]; 

        fread(buffer, sizeof(vec3) * 4, 1, ptr);

        
        memcpy(temp.normal, buffer[0], sizeof(vec3));
        memcpy(temp.vertex1, buffer[1], sizeof(vec3));
        memcpy(temp.vertex2, buffer[2], sizeof(vec3));
        memcpy(temp.vertex3, buffer[3], sizeof(vec3));
        memcpy(temp.colour, (vec3){255, 255, 255}, sizeof(vec3));
        glm_vec3_copy((vec3) {255, 255, 255}, temp.colour);

        scene[1].triangle_list[i] = temp;       

        // Discard attribute 2 byte 
        fseek(ptr, 2, SEEK_CUR);
    }

    

    fclose(ptr);

    // Build save or maybe load bvh?


    srand(time(NULL));

    // Image
    unsigned char* frameData = malloc(OUTPUT_IMAGE_WIDTH * OUTPUT_IMAGE_HEIGHT * 3 * sizeof(char));

    unsigned int index = 0;

    for (int y = OUTPUT_IMAGE_HEIGHT - 1; y >= 0; y--) {
        for (unsigned int x = 0; x < OUTPUT_IMAGE_WIDTH; x++) {
            vec3 pixel_colour = {0, 0, 0};

            for (int i = 0; i < AA_SAMPLES; i++) {
                // Implement BVH here
                double u = ((double) x + randomDouble()) / (OUTPUT_IMAGE_WIDTH - 1);
                double v = ((double) y + randomDouble()) / (OUTPUT_IMAGE_HEIGHT - 1);

                vec3 ray_direction;
                glm_vec3_copy(lower_left_corner, ray_direction);
                vec3 scaled_horizontal, scaled_vertical;
                glm_vec3_scale(horizontal, u, scaled_horizontal);
                glm_vec3_scale(vertical, v, scaled_vertical);
                glm_vec3_add(ray_direction, scaled_horizontal, ray_direction);
                glm_vec3_add(ray_direction, scaled_vertical, ray_direction);
                glm_vec3_sub(ray_direction, origin, ray_direction);

                glm_vec3_normalize(ray_direction);

                vec3 temp;

                trace(origin, ray_direction, &temp, BOUNCE_DEPTH, scene);

                pixel_colour[0] += temp[0];
                pixel_colour[1] += temp[1];
                pixel_colour[2] += temp[2];
                
            }

            glm_vec3_clamp(pixel_colour, 0 , AA_SAMPLES * 255);
            frameData[index] = (unsigned char) (pixel_colour[0] / AA_SAMPLES);
            frameData[index + 1] = (unsigned char) (pixel_colour[1] / AA_SAMPLES);
            frameData[index + 2] = (unsigned char) (pixel_colour[2] / AA_SAMPLES);

            index += 3;
        }
    }

    printf("Completed, saving image...\n");
    stbi_write_png("output.png", OUTPUT_IMAGE_WIDTH, OUTPUT_IMAGE_HEIGHT, 3, frameData, 0);

    free(frameData);

    return 0;
}

void trace(vec3 ray_origin, vec3 ray_direction, vec3* out_pixel_colour, unsigned int depth, object* scene) {
    
    if (depth < 0){
        memcpy(out_pixel_colour, (vec3) {0, 0, 0}, sizeof(vec3));
        return;
    }

    // loop over all the triangles in 1 object
    // issue with overriding
    // look at triangle z index
    vec3 pixel_colour;
    bool changed = false;

    for(int i = 0; i < scene->list_size; i++){
        vec3 vertex1, vertex2, vertex3;

        triangle target_triangle = scene->triangle_list[i];

        glm_vec3_copy(target_triangle.vertex1, vertex1);
        glm_vec3_copy(target_triangle.vertex2, vertex2);
        glm_vec3_copy(target_triangle.vertex3, vertex3);
        
        float distance = -1.0f;
        vec3 intersection;

        ray_triangle_intersection(ray_origin, ray_direction, vertex1, vertex2, vertex3, &distance, &intersection);

        if (distance > 0.0f) {        
            vec3 reflect_ray, temp, surface_normal;

            // True Lambertian 

            random_unit_vector(temp);

            // reading surface normal
            //compute_surface_normal(vertex1, vertex2, vertex3, surface_normal);
            memcpy(surface_normal, target_triangle.normal, sizeof(vec3));

            glm_vec3_add(surface_normal, temp, surface_normal);

            glm_vec3_normalize(surface_normal);

            glm_vec3_reflect(ray_direction, surface_normal, reflect_ray);

            // recursive call
            trace(intersection, reflect_ray, &pixel_colour, depth - 1, scene);
            
            float gamut = 0.5;
            pixel_colour[0] *= gamut;
            pixel_colour[1] *= gamut;
            pixel_colour[2] *= gamut;


            if (changed == false){
                memcpy(out_pixel_colour, pixel_colour, sizeof(vec3));
                changed = true;
            }
            
        }
    }
    if (changed == false){
        float t = 0.5f * (ray_direction[1] + 1.0f);
        pixel_colour[0] = (1.0f - t) * 255.0f + t * 128.0f; // Red component
        pixel_colour[1] = (1.0f - t) * 255.0f + t * 178.0f; // Green component
        pixel_colour[2] = (1.0f - t) * 255.0f + t * 255.0f; // Blue component

        memcpy(out_pixel_colour, pixel_colour, sizeof(vec3));
    }
    

    
}

void compute_surface_normal(vec3 A, vec3 B, vec3 C, vec3 normal) {
    vec3 AB, AC, crossProd;

    // Compute edge vectors AB and AC
    glm_vec3_sub(B, A, AB);
    glm_vec3_sub(C, A, AC);

    // Compute the cross product of AB and AC to get the normal
    glm_vec3_cross(AB, AC, crossProd);

    // Normalize the normal vector
    glm_vec3_normalize(crossProd);

    // Store the result in the provided normal vector
    glm_vec3_copy(crossProd, normal);
}

void random_unit_vector(vec3 dest) {
    // Terrible method
     while (true){
        vec3 temp;
        temp[0] = randomRange(-1,1);
        temp[1] = randomRange(-1,1);
        temp[2] = randomRange(-1,1);

        glm_vec3_normalize(temp);
        if (temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2] >= 1){
            glm_vec3_normalize(temp);
            glm_vec3_copy(temp, dest);
            return ;
        } 
    }
}

double randomRange(int min, int max){
	return min + (max-min)*randomDouble();
}

double randomDouble(){
	return rand() / (RAND_MAX + 1.0);
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

void ray_box_intersection(vec3 ray_origin, vec3 ray_direction, vec3 box_min, vec3 box_max, float* tmin, float* tmax, vec3* out_intersection) {
    const float epsilon = 0.000001;

    vec3 inv_direction;
    glm_vec3_divs(ray_direction, 1.0f, inv_direction); // inverse ray

    vec3 tbot, ttop;
    glm_vec3_mulv(box_min, inv_direction, tbot);
    glm_vec3_mulv(box_max, inv_direction, ttop);

    // Calculate tmin and tmax for each component
    for (int i = 0; i < 3; ++i) {
        tmin[i] = glm_min(tbot[i], ttop[i]);
        tmax[i] = glm_max(tbot[i], ttop[i]);
    }

    *tmin = glm_max(glm_max(tmin[0], tmin[1]), tmin[2]);
    *tmax = glm_min(glm_min(tmax[0], tmax[1]), tmax[2]);

    // Check if intersection is valid
    if (*tmin > *tmax || *tmax < 0) {
        *tmin = -1.0f;
        *tmax = -1.0f;
        return;
    }

    // Calculate intersection point
    glm_vec3_scale(ray_direction, *tmin, *out_intersection);
    glm_vec3_add(ray_origin, *out_intersection, *out_intersection);
}


