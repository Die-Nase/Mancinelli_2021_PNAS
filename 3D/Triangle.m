classdef Triangle < handle
    properties
        id
        verticies
        area
        angles
        l_edges
        v1
        v2
        v3
        f_center
        f_normal
        neighbours
    end
    methods
        function obj = Triangle(id,vertex1,vertex2,vertex3)
            obj.id = id;
            obj.verticies = [vertex1,vertex2,vertex3];
            update(obj);
        end
        
        function update(obj) %vectorize
            for i = 1:length(obj)
                vertex1 = obj(i).verticies(1).xyz;
                vertex2 = obj(i).verticies(2).xyz;
                vertex3 = obj(i).verticies(3).xyz;

                %area
                obj(i).area = 1/2*norm(cross(vertex2-vertex1,vertex3-vertex1));

                %vectors aka edges
                obj(i).v1 = vertex2-vertex1;
                obj(i).v2 = vertex3-vertex2;
                obj(i).v3 = vertex1-vertex3;

                %length of edges
                obj(i).l_edges(1)= norm(obj(i).v1);
                obj(i).l_edges(2)= norm(obj(i).v2);
                obj(i).l_edges(3)= norm(obj(i).v3);

                %triangles
                obj(i).angles(1)=acos(dot(obj(i).v1/obj(i).l_edges(1),-obj(i).v3/obj(i).l_edges(3)));
                obj(i).angles(2)=acos(dot(-obj(i).v1/obj(i).l_edges(1),obj(i).v2/obj(i).l_edges(2)));
                obj(i).angles(3)=pi-(obj(i).angles(1)+obj(i).angles(2));

                %incenter
                a = obj(i).l_edges(2);
                b = obj(i).l_edges(3);
                c = obj(i).l_edges(1);
                obj(i).f_center(1)=(a*vertex1(1)+b*vertex2(1)+c*vertex3(1))/(a+b+c);
                obj(i).f_center(2)=(a*vertex1(2)+b*vertex2(2)+c*vertex3(2))/(a+b+c);
                obj(i).f_center(3)=(a*vertex1(3)+b*vertex2(3)+c*vertex3(3))/(a+b+c);

                %normal vector
                obj(i).f_normal=cross(obj(i).v1,obj(i).v2)/norm(cross(obj(i).v1,obj(i).v2));
            end
        end
        
        function set.neighbours(obj, neighbours)
            obj.neighbours = neighbours;
        end
    end
end