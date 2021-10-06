%% disclaimer
% The curvature calculation in this script uses a modified version of a 
% script published on the matlab fileexchange by Alireza Dastan, 
% which in tutn is based on a paper from Meyer et al. from 2003
%
% citations:
% Alireza Dastan (2021). Gaussian and mean curvatures calculation on a
% triangulated 3d surface (https://www.mathworks.com/matlabcentral/
% fileexchange/61136-gaussian-and-mean-curvatures-calculation-on-a-triangulated-3d-surface),
% MATLAB Central File Exchange. Retrieved October 6, 2021. 
% 
% Meyer, M., Desbrun, M., Schröder, P., & Barr, A. H. (2003). 
% Discrete differential-geometry operators for triangulated 2-manifolds. 
% In Visualization and mathematics III (pp. 35-57). Springer Berlin Heidelberg.
% Available at : http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.24.3427&rep=rep1&type=pdf
%%
classdef Vertex < matlab.mixin.SetGet & matlab.mixin.Copyable
    properties
        id
        
        %% positional properties
        x
        y
        z
        xyz
        theta
        rho
        L
        nnk1
        nnk2
        attached_triangles
        
        %% individual properties
        area
        area0
        mc
        fenergy        
        Ibar_values
        Nbar_values
        Ibar_S
        Nbar_S
        rho_actin
        
        
        
        %% shared properties
        rigidity
        stretch_modulus
        dh
        D
        actin_dh
        Ibar_area
        Nbar_area
        Ibar_c0
        Nbar_c0
        Ibar_mu
        Nbar_mu
        Ibar_rigidity
        Nbar_rigidity
    end
   methods
        function obj = Vertex(id,theta,rho,L)
           obj.id = id;
           obj.theta = theta;
           obj.rho = rho;
           obj.L = L;
           [obj.x,obj.y,obj.z] = pol2cart(theta,rho,L);
        end
        
%         function S_Ibar = get.S_Ibar(obj)
%             S_Ibar = obj.Ibar_area*obj.Ibars/obj.area;
%         end
        function update_protein_saturation(obj,sp)
            for i=1:length(obj)
                obj(i).Ibar_S = obj(i).Ibar_values*sp.Ibar.area/obj(i).area;
                obj(i).Nbar_S = obj(i).Nbar_values*sp.Nbar.area/obj(i).area;
                if obj(i).Ibar_S > 1 || obj(i).Nbar_S > 1
                    error('membrane overflow')
                end
            end
        end

        function xyz = get.xyz(obj)
           xyz = [obj.x,obj.y,obj.z];
        end
        
        function set.nnk1(obj,neighbours)
           obj.nnk1 = neighbours;
        end
        
        function set.nnk2(obj,neighbours)
           obj.nnk2 = neighbours;
        end

        function set.attached_triangles(obj,triangles)
           obj.attached_triangles = triangles;
        end

        function [area,mc] = update(obj)
           mc_vec = [0,0,0];
           n_vec = [0,0,0];
           area = 0;
           point_id = [1,2,3];
           for j = 1:length(obj.attached_triangles)
               T = obj.attached_triangles(j);
               k = point_id(T.verticies == obj);
               % curvature vector
               if k == 1
                   mc_vec= mc_vec+(T.v1/tan(T.angles(3))-T.v3/tan(T.angles(2)));
               elseif k == 2
                   mc_vec= mc_vec+(T.v2/tan(T.angles(1))-T.v1/tan(T.angles(3)));
               elseif k == 3
                   mc_vec= mc_vec+(T.v3/tan(T.angles(2))-T.v2/tan(T.angles(1)));
               end

               %area
                if(T.angles(k)>=pi/2)
                    area = area+T.area/2;
                else
                    if (any(T.angles(:)>=pi/2))
                        area = area+T.area/4;
                    else
                        area2add=0;
                        for m=1:3
                            if m~=k
                                ll=m+1;
                                if ll==4       %% p1==>l2   ,p2==>l3   ,p3==>l1    
                                    ll=1;
                                end
                                area2add=area2add+(T.l_edges(ll)^2/tan(T.angles(m)));
                            end
                        end
                        area=area+area2add/8;
                    end
                end
                wi=1/norm([T.f_center(1)-obj.x,T.f_center(2)-obj.y,T.f_center(3)-obj.z]);
                n_vec=n_vec+wi*T.f_normal;
           end
            mc_vec=mc_vec/4/area;
            n_vec = n_vec/norm(n_vec);
            %sign of mean curvature
            if dot(mc_vec,n_vec) <0
                mc = -norm(mc_vec);
            else
                mc = norm(mc_vec);
            end
            
        end
        
        function move(obj,dh)
                obj.rho = obj.rho + dh;
                [obj.x,obj.y,obj.z] = pol2cart(obj.theta,obj.rho,obj.L);
        end
   end
end