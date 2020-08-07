function [x_remap, y_remap] = electrode_remap(x,y,rois,roi_1_const,roi_midpt,roi_id,c,pa)
    n_remap_steps = 120;
    step_size = 1/n_remap_steps;
    trig_debug = false;
    remap_break_prob = 1;
    
    %remap electrode back into roi
    x_remap = x;
    y_remap = y;

    % Get roi at current point
    diff = (pa.x - x_remap).^2 + (pa.y - y_remap).^2;
    [~,min_idx] = min(diff);
    roi_id(min_idx);
    elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
    
    elec_xy_roi_midpt = roi_midpt(strcmp(rois,roi_1_const),:);
    
    if (~isempty(elec_xy_roi_midpt))
        
        x_remap0 = x_remap;
        y_remap0 = y_remap;
        if (trig_debug)
            plot(x_remap0,y_remap0,'blueo');
            text(x_remap0,y_remap0,roi_1_const,'fontsize',5,'color',0.2*[1 1 1]);
        end

        t_tmp = 0;
        if (~ strcmp(roi_1_const,elec_roi))
            for ir = 1:n_remap_steps
                t_tmp = t_tmp + step_size;
                x_remap = x_remap*(1-t_tmp) + t_tmp * elec_xy_roi_midpt(1);
                y_remap = y_remap*(1-t_tmp) + t_tmp * elec_xy_roi_midpt(2);

                diff = (pa.x - x_remap).^2 + (pa.y - y_remap).^2;
                [~,min_idx] = min(diff);
                roi_id(min_idx);
                elec_roi = c.struct_names{c.table(:,5) == roi_id(min_idx)};
                if strcmp(roi_1_const,elec_roi)
                    %to_break = ( rand([1,1]) <= remap_break_prob );
                    to_break = true;
                    if (to_break)
                        break;
                    end
                end

            end
        end
%         if (~strcmp(elec_roi,'insula'))
%             plot(x_remap,y_remap,'.','Color',...
%                 col_elec_border,'MarkerSize',elec_xy_all_c(ip,6)*border_fac)
%             plot(x_remap,y_remap,'.','Color',...
%                 elec_xy_all_c(ip,3:5),'MarkerSize',elec_xy_all_c(ip,6))
%         end
%         if (trig_debug)
%             plot(x_remap,y_remap,'blue.');
%             plot([x_remap0 x_remap],[y_remap0 y_remap],'blue--')
%             %text(x_remap,y_remap,elec_roi);
%         end
    end

end