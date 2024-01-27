#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <cmath>
#include <limits>
typedef unsigned char RGB[3];

parser::Vec3f cross_product(const parser::Vec3f& edge1, const parser::Vec3f& edge2) {
    parser::Vec3f result;
    result.x = edge1.y * edge2.z - edge1.z * edge2.y;
    result.y = edge1.z * edge2.x - edge1.x * edge2.z;
    result.z = edge1.x * edge2.y - edge1.y * edge2.x;
    return result; // normal of the triangle ?? should I NORMALIZE it?
}

parser::Vec3f multiply_vectors(const parser::Vec3f vec1, const parser::Vec3f vec2){
    parser::Vec3f result;
    result.x = vec1.x * vec2.x;
    result.y = vec1.y * vec2.y;
    result.z = vec1.z * vec2.z;
    return result;
}
parser::Vec3f multiply(const parser::Vec3f& vec, float num){
    parser::Vec3f result;
    result.x = vec.x * num;
    result.y = vec.y * num;
    result.z = vec.z * num;
    return result;
}
parser::Vec3f add(const parser::Vec3f& vec1, const parser::Vec3f& vec2){
    parser::Vec3f result;
    result.x = vec1.x + vec2.x;
    result.y = vec1.y + vec2.y;
    result.z = vec1.z + vec2.z;
    return result;
}

parser::Vec3f subtract(const parser::Vec3f& vec1, const parser::Vec3f& vec2){
    parser::Vec3f result;
    result.x = vec1.x - vec2.x;
    result.y = vec1.y - vec2.y;
    result.z = vec1.z - vec2.z;
    return result;
}
parser::Vec3f divide(const parser::Vec3f& vec, float num){
    parser::Vec3f result = vec;
    if(num != 0){
        result.x = vec.x / num;
        result.y = vec.y / num;
        result.z = vec.z / num;
    }
    return result;

}
parser::Vec3f normalize(const parser::Vec3f& vec){
    parser::Vec3f norm;
    float length = sqrt(vec.x*vec.x + vec.y * vec.y + vec.z*vec.z);
    if (length > 0) return divide(vec,length);
    return vec;
}
float dot_product(const parser::Vec3f& ray_direction, const parser::Vec3f& triangle_normal){
    float result = ray_direction.x * triangle_normal.x + ray_direction.y * triangle_normal.y + ray_direction.z* triangle_normal.z;
    return result; // if sign is positive triangle ignored do not calculate intersection
}
float determinant(float (*arr)[3]){
    //det(A) = a(ei − fh) − b(di − fg) + c(dh − eg)
    float result;
    result = arr[0][0]* ((arr[1][1] * arr[2][2]) - arr[1][2] * arr[2][1]) - 
             arr[0][1]* ((arr[1][0] * arr[2][2]) - arr[1][2] * arr[2][0]) + 
             arr[0][2]* ((arr[1][0] * arr[2][1]) - arr[1][1] * arr[2][0]);
    return result;

}
 


bool triangle_instersection_point(const parser::Vec3f& ray_direction, const parser::Vec3f& position,const parser::Vec3f& a, const parser::Vec3f& b, const parser::Vec3f& c,float& min, bool& inter){
    
    float matrix_A[3][3] = {
        {a.x - c.x, b.x - c.x, -ray_direction.x},
        {a.y-c.y, b.y - c.y, -ray_direction.y},
        {a.z-c.z, b.z - c.z, -ray_direction.z}
    };
    float matrix_alpha[3][3]= {
        {position.x - c.x, b.x - c.x, -ray_direction.x},
        {position.y-c.y, b.y - c.y, -ray_direction.y},
        {position.z-c.z, b.z - c.z, -ray_direction.z}
    };
    float matrix_beta[3][3] = {
        {a.x - c.x, position.x - c.x, -ray_direction.x},
        {a.y-c.y, position.y-c.y, -ray_direction.y},
        {a.z-c.z, position.z-c.z, -ray_direction.z}
    };
    float matrix_tt[3][3] = {
        {a.x - c.x, b.x - c.x, position.x - c.x},
        {a.y-c.y, b.y - c.y, position.y - c.y},
        {a.z-c.z, b.z - c.z, position.z - c.z}
    };
    float alpha,beta,gama,t;
    float det_A = determinant(matrix_A);
    alpha = determinant(matrix_alpha)/ det_A;
    beta = determinant(matrix_beta)/ det_A;
    t = determinant(matrix_tt)/ det_A;
    gama = 1 - (alpha + beta);
    if((alpha >= 0 && alpha <= 1) && (beta >= 0 && beta <= 1) && (gama >= 0 && gama <= 1)){
        if(t > 0 && t < min) {
            min = t;
            inter = true;
            return true;
        } 
    }
    return false;

}

float length(const parser::Vec3f& vec){
    float result = sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
    return result;
}

parser::Vec3i clamp(parser::Vec3f& vec){
    parser::Vec3i result;
    if(vec.x > 255){
        vec.x = 255;
    }
    if(vec.y > 255){
        vec.y = 255;

    } 
    if(vec.z > 255){
        vec.z = 255;
    } 
    result.x = round(vec.x);
    result.y = round(vec.y);
    result.z = round(vec.z);

    return result;
   
}

parser::Vec3f diffuse_shading(const parser::Vec3f& normalized_l, const float& light_dist, const parser::Vec3f& normalized_n, const parser::Material& material, const parser::PointLight& point_light){
    parser::Vec3f diffuse_reflectance, intensity_L, result;
    float cos = std::max((float) 0, dot_product(normalized_n,normalized_l));
    diffuse_reflectance = material.diffuse;
    intensity_L = divide(point_light.intensity, (light_dist*light_dist));
    result = multiply(multiply_vectors(intensity_L, diffuse_reflectance), cos);
    return result;
}

parser::Vec3f specular_shading(const parser::Vec3f& x, parser::Vec3f& normalized_l, const float& light_dist, const parser::Vec3f& normalized_n, const parser::Vec3f& position, const parser::Material& material, const parser::PointLight& point_light){
    parser::Vec3f specular_reflectance, intensity_L, v, h, result;
    float cos;
    v = subtract(position, x);
    h = normalize(add(normalized_l,normalize(v)));
    specular_reflectance = material.specular;
    intensity_L = divide(point_light.intensity, (light_dist*light_dist));
    cos = std::max((float)0, dot_product(normalized_n,h));
    result = multiply(multiply_vectors(intensity_L,specular_reflectance), pow(cos, material.phong_exponent));
    return result;
}

bool is_shadow(const parser::Vec3f& x, const parser::Vec3f& l, const float& light_distance, const parser::Scene& scene, const parser::Vec3f& normalized_n, const parser::Material& material){ 
    float t = std::numeric_limits<float>::infinity(); //???
    bool intersect = false;
    float epsilon = scene.shadow_ray_epsilon;
    float object_distance;
    parser::Vec3f st, object;
    parser::Vec3f intersection_point = add(x, multiply(normalized_n, epsilon));
    parser::Vec3f normalized_l = normalize(l);

    for (int i=0; i < scene.triangles.size(); i++){
        t = std::numeric_limits<float>::infinity();
        intersect = false;
        parser::Triangle triangle = scene.triangles[i];
        int vertexIndex1 = triangle.indices.v0_id - 1; // Vertex indices are 1-based in the XML
        int vertexIndex2 = triangle.indices.v1_id - 1;
        int vertexIndex3 = triangle.indices.v2_id - 1;
        //vertex coordinates
        parser::Vec3f vertex1 = scene.vertex_data[vertexIndex1];
        parser::Vec3f vertex2 = scene.vertex_data[vertexIndex2];
        parser::Vec3f vertex3 = scene.vertex_data[vertexIndex3];
        triangle_instersection_point(l,intersection_point,vertex1,vertex2,vertex3,t,intersect);
        if (intersect && t > 0){
            st = add(intersection_point,multiply(normalized_l,t));
            object = subtract(st,x);
            object_distance = length(object);
            if (object_distance < light_distance){
                return true;
            }
        }
    }

    for(int k = 0; k < scene.spheres.size(); k++){
        int center_id;
        parser::Sphere sphere = scene.spheres[k];
        float radius;
        parser::Vec3f center_position;
        float a,b,c,dis,t1,t2,t;
        center_id = sphere.center_vertex_id - 1;
        radius = sphere.radius;
        center_position = scene.vertex_data[center_id];
        parser::Vec3f sub_inter_center = subtract(intersection_point,center_position);
        a = dot_product(l,l);
        b = dot_product(multiply(l,2),sub_inter_center);
        c = dot_product(sub_inter_center,sub_inter_center) - radius*radius;
        dis = b*b - 4*a*c;
        if(dis >= 0){
            t1 = (-b - sqrt(dis))/(2*a); // ??? diğer t
            t2 = (-b + sqrt(dis))/(2*a);
            if (t1 > 0) {
                st = add(intersection_point,multiply(normalized_l,t1));
                object = subtract(st,x);
                object_distance = length(object);
                if (object_distance < light_distance){
                    return true;
                }
            }
            
            if (t2 > 0){
                st = add(intersection_point,multiply(normalized_l,t2));
                object = subtract(st,x);
                object_distance = length(object);
                if (object_distance < light_distance){
                    return true;
                }
            }
        }
    }

    for(int k = 0; k <  scene.meshes.size(); k++){
        parser::Mesh mesh = scene.meshes[k];
        for(int n = 0; n < mesh.faces.size(); n++){
            t = std::numeric_limits<float>::infinity();
            intersect = false;
            parser::Vec3f triangle_intersection;
            int vertexIndex1 = mesh.faces[n].v0_id - 1; // Vertex indices are 1-based in the XML
            int vertexIndex2 = mesh.faces[n].v1_id - 1;
            int vertexIndex3 = mesh.faces[n].v2_id - 1;
            //vertex coordinates
            parser::Vec3f vertex1 = scene.vertex_data[vertexIndex1];
            parser::Vec3f vertex2 = scene.vertex_data[vertexIndex2];
            parser::Vec3f vertex3 = scene.vertex_data[vertexIndex3];
            triangle_instersection_point(l,intersection_point,vertex1,vertex2,vertex3,t,intersect);
            if (intersect && t > 0){
                st = add(intersection_point,multiply(normalized_l,t));
                object = subtract(st,x);
                object_distance = length(object);
                if (object_distance < light_distance){
                    return true;
                }
            }
        }
    }

    return false;
}

parser::Vec3f compute_color(const parser::Scene& scene, const parser::Vec3f& ray_direction, const parser::Vec3f& position, int& depth);
parser::Vec3f apply_shading(const parser::Vec3f& x, const parser::Vec3f ray_direction, const parser::Scene& scene, const parser::Vec3f& position, const parser::Vec3f& normalized_n, const parser::Material& material, int& depth){ 
    parser::Vec3f ambient_color, diffuse_color, specular_color, l, normalized_l;
    parser::Vec3f ambient_reflectance = material.ambient;
    parser::Vec3f ambient_radience = scene.ambient_light;
    float ligth_distance;
    ambient_color = multiply_vectors(ambient_reflectance,ambient_radience);

    if (material.is_mirror){ //mirror
        depth++;
        parser::Vec3f w_r, w_o, pos;
        w_o = normalize(multiply(ray_direction, -1));
        w_r = normalize(add(multiply(w_o, -1), multiply(normalized_n,dot_product(normalized_n,w_o)*2)));
        pos = add(x, multiply(normalized_n, scene.shadow_ray_epsilon));
        ambient_color = add(ambient_color, multiply_vectors(compute_color(scene, w_r, pos, depth), material.mirror));
    }

    for (int i=0; i < scene.point_lights.size(); i++){
        parser::PointLight point_light = scene.point_lights[i];
        l = subtract(point_light.position, x);
        ligth_distance = length(l);
        normalized_l = normalize(l);
        if (!is_shadow(x,l,ligth_distance,scene,normalized_n,material)){
            diffuse_color = diffuse_shading(normalized_l,ligth_distance,normalized_n,material,point_light);
            specular_color = specular_shading(x,normalized_l, ligth_distance, normalized_n, position, material, point_light);
            ambient_color = add(ambient_color, diffuse_color);
            ambient_color = add(ambient_color, specular_color);
        }
    }

    return ambient_color;

}


parser::Vec3f compute_color(const parser::Scene& scene, const parser::Vec3f& ray_direction, const parser::Vec3f& position, int& depth){ 
    parser::Vec3f pixel_color;
    if (depth > scene.max_recursion_depth){
        pixel_color.x = 0;
        pixel_color.y = 0;
        pixel_color.z = 0;
        return pixel_color;
    }

    float t_min = std::numeric_limits<float>::infinity(); //???
    bool intersect = false;
    parser::Vec3f normal;
    int material_id;
    for(int k = 0; k < scene.triangles.size(); k++){
        parser::Triangle triangle = scene.triangles[k];
        int vertexIndex1 = triangle.indices.v0_id - 1; // Vertex indices are 1-based in the XML
        int vertexIndex2 = triangle.indices.v1_id - 1;
        int vertexIndex3 = triangle.indices.v2_id - 1;
        //vertex coordinates
        parser::Vec3f vertex1 = scene.vertex_data[vertexIndex1];
        parser::Vec3f vertex2 = scene.vertex_data[vertexIndex2];
        parser::Vec3f vertex3 = scene.vertex_data[vertexIndex3];
        parser::Vec3f edge1 = subtract(vertex2,vertex1);
        parser::Vec3f edge2 = subtract(vertex3,vertex1);
        parser::Vec3f triangle_normal = normalize(cross_product(edge1, edge2));
        float back_face = dot_product(ray_direction,triangle_normal);
        if(back_face < 0){// Back_face_culling
            bool tri = triangle_instersection_point(ray_direction,position,vertex1,vertex2,vertex3,t_min,intersect);
            if(tri == true){
                material_id = triangle.material_id;
                normal = triangle_normal;
            }
        }
    }

    for(int k = 0; k < scene.spheres.size(); k++){
        int center_id;
        parser::Sphere sphere = scene.spheres[k];
        float radius;
        parser::Vec3f center_position;
        float a,b,c,t1,t2,dis;
        center_id = sphere.center_vertex_id - 1;
        radius = sphere.radius;
        center_position = scene.vertex_data[center_id];
        parser::Vec3f sub_pos_center = subtract(position,center_position);
        a = dot_product(ray_direction,ray_direction);
        b = dot_product(multiply(ray_direction,2),sub_pos_center);
        c = dot_product(sub_pos_center,sub_pos_center) - radius*radius;
        dis = b*b - 4*a*c;
        if(dis >= 0){
            t1 = (-b - sqrt(dis))/(2*a); // ??? diğer t
            if(t1 < t_min && t1 > 0){
                t_min = t1;
                intersect = true;
                material_id = sphere.material_id;
                normal = divide(subtract(add(position, multiply(ray_direction,t_min)), center_position), radius);
            }
        } 
    }

    for(int k = 0; k <  scene.meshes.size(); k++){
        parser::Mesh mesh = scene.meshes[k];
        for(int n = 0; n < mesh.faces.size(); n++){
            parser::Vec3f triangle_intersection;
            int vertexIndex1 = mesh.faces[n].v0_id - 1; // Vertex indices are 1-based in the XML
            int vertexIndex2 = mesh.faces[n].v1_id - 1;
            int vertexIndex3 = mesh.faces[n].v2_id - 1;
            //vertex coordinates
            parser::Vec3f vertex1 = scene.vertex_data[vertexIndex1];
            parser::Vec3f vertex2 = scene.vertex_data[vertexIndex2];
            parser::Vec3f vertex3 = scene.vertex_data[vertexIndex3];
            parser::Vec3f edge1 = subtract(vertex2,vertex1);
            parser::Vec3f edge2 = subtract(vertex3,vertex1);
            parser::Vec3f triangle_normal = normalize(cross_product(edge1, edge2));
            float back_face = dot_product(ray_direction,triangle_normal);
            if (back_face < 0){
                bool tri = triangle_instersection_point(ray_direction,position,vertex1,vertex2,vertex3,t_min,intersect);
                if(tri == true){
                    material_id = mesh.material_id;
                    normal = triangle_normal;
                }
            }
        }
    }

    if(intersect == true){
        parser::Vec3f intersection_point = add(position, multiply(ray_direction,t_min));
        parser::Material material = scene.materials[material_id-1];
        parser::Vec3f normalized_n = normalize(normal);
        return apply_shading(intersection_point, ray_direction, scene, position, normalized_n, material, depth);
    }

    else if (depth == 0){
        pixel_color.x = scene.background_color.x;
        pixel_color.y = scene.background_color.y;
        pixel_color.z = scene.background_color.z;
        return pixel_color;
    }
    else {
        pixel_color.x = 0;
        pixel_color.y = 0;
        pixel_color.z = 0;
        return pixel_color;
    }
    

}




int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;
    scene.loadFromXml(argv[1]);
    for(int a = 0; a < scene.cameras.size(); a++){//??
        parser::Vec3f ray_direction;
        parser::Camera camera = scene.cameras[a];
        float left = camera.near_plane.x;
        float right = camera.near_plane.y;
        float bottom = camera.near_plane.z;
        float top = camera.near_plane.w;//????
        float pixel_width = (right - left)/camera.image_width;
        float pixel_height = (top - bottom)/camera.image_height;
        parser::Vec3f u = cross_product(camera.up, multiply(camera.gaze,-1));
        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];
        int pixel_counter = 0;
        parser::Vec3f m, q, s, pixel_color;
        float s_u, s_v;
        m = add(camera.position, multiply(camera.gaze,camera.near_distance));
        q = add(add(m, multiply(u,left)),multiply(camera.up, top));

        for(int i = 0; i < camera.image_height; i++){ 
            for(int j = 0; j < camera.image_width; j++ ){ //!!!!
                int depth = 0;
                s_u = (j + 0.5)*pixel_width;
                s_v = (i + 0.5)*pixel_height;
                s = add(q, subtract(multiply(u, s_u),multiply(camera.up,s_v)));
                ray_direction = subtract(s,camera.position);
                
                pixel_color = compute_color(scene, ray_direction, camera.position, depth);
                parser::Vec3i color = clamp(pixel_color);

                image[pixel_counter++] = color.x;
                image[pixel_counter++] = color.y;
                image[pixel_counter++] = color.z;

            }
        }
        const char* name = camera.image_name.c_str();
        write_ppm(name, image, camera.image_width, camera.image_height);

    }


}