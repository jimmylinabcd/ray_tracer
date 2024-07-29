#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <time.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include <cglm/cglm.h>
#include <cglm/call.h>
#include "tinycthread.h"
#include <stddef.h>

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define OUTPUT_IMAGE_HEIGHT 1080
#define OUTPUT_IMAGE_WIDTH 1080
#define AA_SAMPLES 2
#define BOUNCE_DEPTH 5
#define NUM_THREADS 4

typedef vec3 mat3[3];
#ifdef __AVX__
typedef CGLM_ALIGN_IF(32) vec4 mat4[4];
#else
typedef CGLM_ALIGN_IF(16) vec4 mat4[4];
#endif

typedef struct {
    vec3 min;
    vec3 max;
} bounding_box;

typedef struct BVH_node{
    bounding_box bdbox;
    struct BVH_node* left;
    struct BVH_node* right;
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

typedef struct {
    unsigned char* frameData;
    int thread_id;
    object* scene;
    vec3 lower_left_corner;
    vec3 horizontal;
    vec3 vertical;
    vec3 origin;
} thread_data;

void ray_triangle_intersection(vec3 ray_origin, vec3 ray_vector, vec3 vertex1, vec3 vertex2, vec3 vertex3, float* t, vec3* out_intersection);
void ray_box_intersection(vec3 ray_origin, vec3 ray_direction, vec3 box_min, vec3 box_max, float* tmin, float* tmax, vec3* out_intersection);
void trace(vec3 ray_origin, vec3 ray_direction, vec3* out_pixel_colour, unsigned int depth, object* scene);
void compute_surface_normal(vec3 A, vec3 B, vec3 C, vec3 normal);
double randomRange(double min, double max);
double randomDouble();
void create_triangle(triangle *tri, vec3 v1, vec3 v2, vec3 v3, vec3 color);
void scale_and_center_model(object* obj, vec3 center_position);
void random_unit_vector(vec3 dest);
int render_thread(void* arg);

int main() {
    printf("Running Ray Tracer\n");

    srand(time(NULL));

    // Scene Description

    object scene[2];

    vec3 origin = {0.0f, 0.0f, 0.0f};
    double focal_length = 1;
    double viewport_height = 2.0;
    double viewport_width = ((double)OUTPUT_IMAGE_WIDTH / (double)OUTPUT_IMAGE_HEIGHT) * viewport_height;

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

    // Floor
    triangle floor_triangle1;
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -1.0f}, floor_triangle1.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -1.0f}, floor_triangle1.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, floor_triangle1.vertex3);
    vec3 temp_normal1;
    compute_surface_normal(floor_triangle1.vertex1, floor_triangle1.vertex2, floor_triangle1.vertex3, temp_normal1);

    triangle floor_triangle2;
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -1.0f}, floor_triangle2.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, floor_triangle2.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, floor_triangle2.vertex3);
    vec3 temp_normal2;
    compute_surface_normal(floor_triangle2.vertex1, floor_triangle2.vertex2, floor_triangle2.vertex3, temp_normal2);

    memcpy(floor_triangle1.normal, temp_normal1, sizeof(vec3));
    memcpy(floor_triangle2.normal, temp_normal2, sizeof(vec3));

    glm_vec3_copy((vec3){255, 255, 255}, floor_triangle1.colour);
    glm_vec3_copy((vec3){255, 255, 255}, floor_triangle2.colour);

    // Ceiling
    triangle ceiling_triangle1;
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -1.0f}, ceiling_triangle1.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -1.0f}, ceiling_triangle1.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -3.0f}, ceiling_triangle1.vertex3);
    vec3 temp_normal3;
    compute_surface_normal(ceiling_triangle1.vertex1, ceiling_triangle1.vertex2, ceiling_triangle1.vertex3, temp_normal3);

    triangle ceiling_triangle2;
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -1.0f}, ceiling_triangle2.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -3.0f}, ceiling_triangle2.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -3.0f}, ceiling_triangle2.vertex3);
    vec3 temp_normal4;
    compute_surface_normal(ceiling_triangle2.vertex1, ceiling_triangle2.vertex2, ceiling_triangle2.vertex3, temp_normal4);

    memcpy(ceiling_triangle1.normal, temp_normal3, sizeof(vec3));
    memcpy(ceiling_triangle2.normal, temp_normal4, sizeof(vec3));

    glm_vec3_copy((vec3){255, 255, 255}, ceiling_triangle1.colour);
    glm_vec3_copy((vec3){255, 255, 255}, ceiling_triangle2.colour);

    // Left wall
    triangle left_wall_triangle1;
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -1.0f}, left_wall_triangle1.vertex1);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -1.0f}, left_wall_triangle1.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, left_wall_triangle1.vertex3);
    vec3 temp_normal5;
    compute_surface_normal(left_wall_triangle1.vertex1, left_wall_triangle1.vertex2, left_wall_triangle1.vertex3, temp_normal5);

    triangle left_wall_triangle2;
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -1.0f}, left_wall_triangle2.vertex1);
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, left_wall_triangle2.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -3.0f}, left_wall_triangle2.vertex3);
    vec3 temp_normal6;
    compute_surface_normal(left_wall_triangle2.vertex1, left_wall_triangle2.vertex2, left_wall_triangle2.vertex3, temp_normal6);

    memcpy(left_wall_triangle1.normal, temp_normal5, sizeof(vec3));
    memcpy(left_wall_triangle2.normal, temp_normal6, sizeof(vec3));

    glm_vec3_copy((vec3){255, 0, 0}, left_wall_triangle1.colour);
    glm_vec3_copy((vec3){255, 0, 0}, left_wall_triangle2.colour);

    // Right wall
    triangle right_wall_triangle1;
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -1.0f}, right_wall_triangle1.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -1.0f}, right_wall_triangle1.vertex2);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, right_wall_triangle1.vertex3);
    vec3 temp_normal7;
    compute_surface_normal(right_wall_triangle1.vertex1, right_wall_triangle1.vertex2, right_wall_triangle1.vertex3, temp_normal7);

    triangle right_wall_triangle2;
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -1.0f}, right_wall_triangle2.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, right_wall_triangle2.vertex2);
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -3.0f}, right_wall_triangle2.vertex3);
    vec3 temp_normal8;
    compute_surface_normal(right_wall_triangle2.vertex1, right_wall_triangle2.vertex2, right_wall_triangle2.vertex3, temp_normal8);

    memcpy(right_wall_triangle1.normal, temp_normal7, sizeof(vec3));
    memcpy(right_wall_triangle2.normal, temp_normal8, sizeof(vec3));

    glm_vec3_copy((vec3){0, 0, 255}, right_wall_triangle1.colour);
    glm_vec3_copy((vec3){0, 0, 255}, right_wall_triangle2.colour);

    // Back wall
    triangle back_wall_triangle1;
    glm_vec3_copy((vec3){-viewport_width / 2, -viewport_height / 2, -3.0f}, back_wall_triangle1.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, back_wall_triangle1.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -3.0f}, back_wall_triangle1.vertex3);
    vec3 temp_normal9;
    compute_surface_normal(back_wall_triangle1.vertex1, back_wall_triangle1.vertex2, back_wall_triangle1.vertex3, temp_normal9);

    triangle back_wall_triangle2;
    glm_vec3_copy((vec3){viewport_width / 2, -viewport_height / 2, -3.0f}, back_wall_triangle2.vertex1);
    glm_vec3_copy((vec3){viewport_width / 2, viewport_height / 2, -3.0f}, back_wall_triangle2.vertex2);
    glm_vec3_copy((vec3){-viewport_width / 2, viewport_height / 2, -3.0f}, back_wall_triangle2.vertex3);
    vec3 temp_normal10;
    compute_surface_normal(back_wall_triangle2.vertex1, back_wall_triangle2.vertex2, back_wall_triangle2.vertex3, temp_normal10);

    memcpy(back_wall_triangle1.normal, temp_normal9, sizeof(vec3));
    memcpy(back_wall_triangle2.normal, temp_normal10, sizeof(vec3));

    glm_vec3_copy((vec3){0, 255, 0}, back_wall_triangle1.colour);
    glm_vec3_copy((vec3){0, 255, 0}, back_wall_triangle2.colour);


    triangle light_triangle1;
    vec3 light_center = {0, (viewport_height / 2) - 0.5, -2.0f};
    float light_size = 1.0f;

    glm_vec3_copy((vec3){light_center[0] - light_size / 2, light_center[1] + light_size / 2, light_center[2]}, light_triangle1.vertex1);
    glm_vec3_copy((vec3){light_center[0] + light_size / 2, light_center[1] + light_size / 2, light_center[2]}, light_triangle1.vertex2);
    glm_vec3_copy((vec3){light_center[0] - light_size / 2, light_center[1] - light_size / 2, light_center[2]}, light_triangle1.vertex3);
    vec3 temp_normal_light1;
    compute_surface_normal(light_triangle1.vertex1, light_triangle1.vertex2, light_triangle1.vertex3, temp_normal_light1);

    triangle light_triangle2;
    glm_vec3_copy((vec3){light_center[0] + light_size / 2, light_center[1] + light_size / 2, light_center[2]}, light_triangle2.vertex1);
    glm_vec3_copy((vec3){light_center[0] - light_size / 2, light_center[1] - light_size / 2, light_center[2]}, light_triangle2.vertex2);
    glm_vec3_copy((vec3){light_center[0] + light_size / 2, light_center[1] - light_size / 2, light_center[2]}, light_triangle2.vertex3);
    vec3 temp_normal_light2;
    compute_surface_normal(light_triangle2.vertex1, light_triangle2.vertex2, light_triangle2.vertex3, temp_normal_light2);

    memcpy(light_triangle1.normal, temp_normal_light1, sizeof(vec3));
    memcpy(light_triangle2.normal, temp_normal_light2, sizeof(vec3));

    // Set the light color (e.g., white)
    glm_vec3_copy((vec3){255, 255, 255}, light_triangle1.colour);
    glm_vec3_copy((vec3){255, 255, 255}, light_triangle2.colour);

    // create test object
    scene[0].id = 0;
    //scene[0].triangle_list = &test_triangle;
    scene[0].triangle_list = malloc(sizeof(triangle) * 12);
    scene[0].triangle_list[0] = floor_triangle1;
    scene[0].triangle_list[1] = floor_triangle2;
    scene[0].triangle_list[2] = ceiling_triangle1;
    scene[0].triangle_list[3] = ceiling_triangle2;
    scene[0].triangle_list[4] = left_wall_triangle1;
    scene[0].triangle_list[5] = left_wall_triangle2;
    scene[0].triangle_list[6] = right_wall_triangle1;
    scene[0].triangle_list[7] = right_wall_triangle2;
    scene[0].triangle_list[8] = back_wall_triangle1;
    scene[0].triangle_list[9] = back_wall_triangle2;
    scene[0].triangle_list[10] = light_triangle1;
    scene[0].triangle_list[11] = light_triangle2;
    scene[0].list_size = 12;
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


    // Load vertecies into memeory
    scene[1].id = 1;
    scene[1].triangle_list = malloc(sizeof(triangle) * number_triangles);
    scene[1].list_size = number_triangles;
    scene[1].root = test_bvh;

    vec3 *verticies = (vec3  *)malloc(number_triangles  * 3 * sizeof(vec3));

	if (verticies == NULL) {
        printf("Memory allocation failed\n");
        fclose(ptr);
        return 1;
    }

    for (unsigned int i = 0; i < number_triangles; i++) {
        // Discard normal vector 12 bytes
        fseek(ptr, 12, SEEK_CUR);

        // Read verticies 3 * 12 bytes
        fread(&verticies[i * 3], sizeof(vec3), 3, ptr);

	    // Discard attribute byte count 2 bytes
        fseek(ptr, 2, SEEK_CUR);
    }

    fclose(ptr);

	printf("Loaded %d triangles\n", number_triangles);

	for (unsigned int i = 0; i < number_triangles; i++) {
        triangle *tri = &scene[1].triangle_list[i];
        glm_vec3_copy(verticies[i * 3], tri->vertex1);
        glm_vec3_copy(verticies[i * 3 + 1], tri->vertex2);
        glm_vec3_copy(verticies[i * 3 + 2], tri->vertex3);

        // Compute normal
        compute_surface_normal(tri->vertex1, tri->vertex2, tri->vertex3, tri->normal);

        // Set color to white
        glm_vec3_copy((vec3){255, 255, 255}, tri->colour);
    }

    //for (unsigned int i = 0; i < number_triangles; i++) {
        //printf("(%f, %f, %f)\n", verticies[i][0], verticies[i][1], verticies[i][2]);
        //printf("(%f, %f, %f)\n", scene[1].triangle_list[i].vertex1[0], scene[1].triangle_list[i].vertex1[1], scene[1].triangle_list[i].vertex1[2]);
        //printf("(%f, %f, %f)\n", scene[1].triangle_list[i].vertex2[0], scene[1].triangle_list[i].vertex2[1], scene[1].triangle_list[i].vertex2[2]);
        //printf("(%f, %f, %f)\n", scene[1].triangle_list[i].vertex3[0], scene[1].triangle_list[i].vertex3[1], scene[1].triangle_list[i].vertex3[2]);

    //}


    vec3 center_position = {0.0f, 0.0f, -0.1f};  // Example center position
    scale_and_center_model(&scene[1], center_position);  // Assuming the knight model is in scene[1]

    // Image
    unsigned char* frameData = malloc(OUTPUT_IMAGE_WIDTH * OUTPUT_IMAGE_HEIGHT * 3 * sizeof(char));

   thread_data t_data[NUM_THREADS];
    thrd_t threads[NUM_THREADS];

    // Divide work among threads
    for (int i = 0; i < NUM_THREADS; i++) {
        t_data[i].frameData = frameData;
        t_data[i].thread_id = i;
        t_data[i].scene = scene;
        glm_vec3_copy(lower_left_corner, t_data[i].lower_left_corner);
        glm_vec3_copy(horizontal, t_data[i].horizontal);
        glm_vec3_copy(vertical, t_data[i].vertical);
        glm_vec3_copy(origin, t_data[i].origin);

        thrd_create(&threads[i], render_thread, &t_data[i]);
    }

    // Wait for all threads to complete
    for (int i = 0; i < NUM_THREADS; i++) {
        thrd_join(threads[i], NULL);
    }

    printf("Completed, saving image...\n");
    stbi_write_png("output.png", OUTPUT_IMAGE_WIDTH, OUTPUT_IMAGE_HEIGHT, 3, frameData, 0);

    free(frameData);

    return 0;
}

void trace(vec3 ray_origin, vec3 ray_direction, vec3* out_pixel_colour, unsigned int depth, object* scene) {
    
    if (depth == 0){
        memcpy(out_pixel_colour, (vec3) {0, 0, 0}, sizeof(vec3));
        return;
    }

    vec3 pixel_colour;
    bool changed = false;

    for (unsigned int obj_idx = 0; obj_idx < 2; obj_idx++) {
        object* obj = &scene[obj_idx];

        for(int i = 0; i < obj->list_size; i++){
            vec3 vertex1, vertex2, vertex3;

            triangle target_triangle = obj->triangle_list[i];

            glm_vec3_copy(target_triangle.vertex1, vertex1);
            glm_vec3_copy(target_triangle.vertex2, vertex2);
            glm_vec3_copy(target_triangle.vertex3, vertex3);
                        
            float distance = -1.0f;
            vec3 intersection;

            ray_triangle_intersection(ray_origin, ray_direction, vertex1, vertex2, vertex3, &distance, &intersection);

            if (distance > 0.0f && distance < FLT_MAX) {        
                vec3 reflect_ray, temp, surface_normal;

                random_unit_vector(temp);

                glm_vec3_copy(target_triangle.normal, surface_normal);

                glm_vec3_add(surface_normal, temp, surface_normal);

                glm_vec3_normalize(surface_normal);

                glm_vec3_reflect(ray_direction, surface_normal, reflect_ray);

                // recursive call
                trace(intersection, reflect_ray, &pixel_colour, depth - 1, scene);

                pixel_colour[0] *= target_triangle.colour[0]/255.0;
                pixel_colour[1] *= target_triangle.colour[1]/255.0;
                pixel_colour[2] *= target_triangle.colour[2]/255.0;
                
                float gamut = 0.6;
                pixel_colour[0] *= gamut;
                pixel_colour[1] *= gamut;
                pixel_colour[2] *= gamut;

                if (changed == false){
                    memcpy(out_pixel_colour, pixel_colour, sizeof(vec3));
                    changed = true;
                }
                
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
        temp[0] = randomRange(-1.0,1.0);
        temp[1] = randomRange(-1.0,1.0);
        temp[2] = randomRange(-1.0,1.0);

        if ((temp[0]*temp[0] + temp[1]*temp[1] + temp[2]*temp[2]) < 1){
            glm_vec3_normalize(temp);
            glm_vec3_copy(temp, dest);
            return ;
        } 
    }
}

double randomRange(double min, double max){
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

void scale_and_center_model(object* obj, vec3 center_position) {
    // Find the bounding box of the object
    vec3 min_bound = {FLT_MAX, FLT_MAX, FLT_MAX};
    vec3 max_bound = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    for (unsigned int i = 0; i < obj->list_size; i++) {
        triangle* tri = &obj->triangle_list[i];
        vec3* vertices[] = {&tri->vertex1, &tri->vertex2, &tri->vertex3};
        for (int j = 0; j < 3; j++) {
            glm_vec3_minv(*vertices[j], min_bound, min_bound);
            glm_vec3_maxv(*vertices[j], max_bound, max_bound);
        }
    }

    // Calculate the size in each dimension
    vec3 size;
    glm_vec3_sub(max_bound, min_bound, size);

    // Find the largest dimension
    float max_size = fmaxf(fmaxf(size[0], size[1]), size[2]);

    float scale_factor = 0.5f / max_size;

    vec3 box_center;
    glm_vec3_add(min_bound, max_bound, box_center);
    glm_vec3_scale(box_center, 0.5f, box_center);

    glm_vec3_scale(box_center, scale_factor, box_center);

    vec3 translation;
    glm_vec3_sub(center_position, box_center, translation);

    // Apply the scaling and translation to each vertex
    for (unsigned int i = 0; i < obj->list_size; i++) {
        triangle* tri = &obj->triangle_list[i];
        vec3* vertices[] = {&tri->vertex1, &tri->vertex2, &tri->vertex3};
        for (int j = 0; j < 3; j++) {
            glm_vec3_scale(*vertices[j], scale_factor, *vertices[j]);
            glm_vec3_add(*vertices[j], translation, *vertices[j]);
        }
    }

    // Recalculate the bounding box after transformation
    vec3 new_min_bound = {FLT_MAX, FLT_MAX, FLT_MAX};
    vec3 new_max_bound = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

    for (unsigned int i = 0; i < obj->list_size; i++) {
        triangle* tri = &obj->triangle_list[i];
        vec3* vertices[] = {&tri->vertex1, &tri->vertex2, &tri->vertex3};
        for (int j = 0; j < 3; j++) {
            glm_vec3_minv(*vertices[j], new_min_bound, new_min_bound);
            glm_vec3_maxv(*vertices[j], new_max_bound, new_max_bound);
        }
    }

    // Print the new bounding box
    printf("New bounding box:\n");
    printf("MIN: %f, %f, %f\n", new_min_bound[0], new_min_bound[1], new_min_bound[2]);
    printf("MAX: %f, %f, %f\n", new_max_bound[0], new_max_bound[1], new_max_bound[2]);
}


void create_triangle(triangle *tri, vec3 v1, vec3 v2, vec3 v3, vec3 color) {
    glm_vec3_copy(v1, tri->vertex1);
    glm_vec3_copy(v2, tri->vertex2);
    glm_vec3_copy(v3, tri->vertex3);
    compute_surface_normal(v1, v2, v3, tri->normal);
    glm_vec3_copy(color, tri->colour);
}

int render_thread(void* arg) {
    thread_data* tdata = (thread_data*)arg;

    int start = (OUTPUT_IMAGE_HEIGHT / NUM_THREADS) * tdata->thread_id;
    int end = (OUTPUT_IMAGE_HEIGHT / NUM_THREADS) * (tdata->thread_id + 1);

    for (int y = end - 1; y >= start; y--) {
        for (unsigned int x = 0; x < OUTPUT_IMAGE_WIDTH; x++) {
            vec3 pixel_colour = {0, 0, 0};

            for (int i = 0; i < AA_SAMPLES; i++) {
                double u = ((double)x + randomDouble()) / (OUTPUT_IMAGE_WIDTH - 1);
                double v = ((double)y + randomDouble()) / (OUTPUT_IMAGE_HEIGHT - 1);

                vec3 ray_direction;
                glm_vec3_copy(tdata->lower_left_corner, ray_direction);
                vec3 scaled_horizontal, scaled_vertical;
                glm_vec3_scale(tdata->horizontal, u, scaled_horizontal);
                glm_vec3_scale(tdata->vertical, v, scaled_vertical);
                glm_vec3_add(ray_direction, scaled_horizontal, ray_direction);
                glm_vec3_add(ray_direction, scaled_vertical, ray_direction);
                glm_vec3_sub(ray_direction, tdata->origin, ray_direction);

                glm_vec3_normalize(ray_direction);

                vec3 temp;
                trace(tdata->origin, ray_direction, &temp, BOUNCE_DEPTH, tdata->scene);

                pixel_colour[0] += temp[0];
                pixel_colour[1] += temp[1];
                pixel_colour[2] += temp[2];
            }

            glm_vec3_clamp(pixel_colour, 0, AA_SAMPLES * 255);
            int index = (y * OUTPUT_IMAGE_WIDTH + x) * 3;
            tdata->frameData[index] = (unsigned char)(pixel_colour[0] / AA_SAMPLES);
            tdata->frameData[index + 1] = (unsigned char)(pixel_colour[1] / AA_SAMPLES);
            tdata->frameData[index + 2] = (unsigned char)(pixel_colour[2] / AA_SAMPLES);
        }
    }

    return 0;
}