close all;
clear;

SubjectsSleep = {'m00019','m00023','m00024','m00026','m00030','m00032','m00035',...
    'm00039','m00043','m00044','m00049','m00079','m00083','m00084'};

Subjects24h = {'m00030','m00043','m00083'};

mkdir('./figures');
mkdir('./figures/T21');

Btxt = {'sleep','24h'};
%Btxt = {'sleep'};

trig_plot_coh_only = true;

for metrici = [1] %[1 5]
    for ib = 1:length(Btxt)
        
        fprintf('[*] metric %i, %s\n',metrici, Btxt{ib});
    
        if (strcmp(Btxt{ib},'24h'))
            T30 = load('~/Dropbox/24hours/behavior/m00030_Annot.mat');
            T43 = load('~/Dropbox/24hours/behavior/m00043_Annot.mat');
            T83 = load('~/Dropbox/24hours/behavior/m00083_Annot.mat');
            Behs = T83.AnnotLabels;
            IB = [1 2 4 10];
%             [T30.AnnotLabels,T43.AnnotLabels,T83.AnnotLabels]
% 
%             ans =
% 
%               10×3 cell array
% 
%                 {'ARM MOVEMENT'      }    {'BODY MOVEMENT'     }    {'BODY MOVEMENT'     }
%                 {'CONTACT'           }    {'CONTACT'           }    {'CONTACT'           }
%                 {'EATING'            }    {'EATING'            }    {'EATING'            }
%                 {'HEAD MOVEMENT'     }    {'HEAD MOVEMENT'     }    {'HEAD MOVEMENT'     }
%                 {'LEG MOVEMENT'      }    {'PATIENT IS TALKING'}    {'PATIENT IS TALKING'}
%                 {'PATIENT IS TALKING'}    {'QUIET'             }    {'QUIET'             }
%                 {'SLEEP'             }    {'SLEEP'             }    {'SLEEP'             }
%                 {'SOMEONE IS TALKING'}    {'SOMEONE IS TALKING'}    {'SOMEONE IS TALKING'}
%                 {'VIDEO GAMES'       }    {'VIDEO GAMES'       }    {'VIDEO GAMES'       }
%                 {'WATCH TV'          }    {'WATCH TV'          }    {'WATCH TV'          }
            %T.AnnotLabels
        else
            Behs = {NaN};
            IB = 1:length(Behs);
        end
        

        for ib2 = IB
            
            Beh_txt = Behs{ib2};
        
            % Sleep subjects
            if (any(isnan(Behs{ib2})))
                state_24h = false;
                infname = sprintf('./cache/plot_%s_sw_train_metric-%i.mat',Btxt{ib},metrici);
                infname_perm = sprintf('./cache/plot_%s_sw_train_metric-%i_perm.mat',Btxt{ib},metrici);
                pos = [0 0 800 200];
                n_subj = length(SubjectsSleep);
            else
                state_24h = true;
                btxt = replace(Behs{ib2},' ','_');
                pos = [0 0 324 200];
                infname = sprintf('./cache/plot_%s_sw_train_metric-%i_%s.mat',Btxt{ib},metrici,btxt);
                infname_perm = sprintf('./cache/plot_%s_sw_train_metric-%i_%s_perm.mat',Btxt{ib},metrici,btxt);
                n_subj = length(Subjects24h);
            end
            
            Ca = load(infname);
            Ca_perm = load(infname_perm);
            
            % Accuracy threshold
            n_sub_perm = length(Ca_perm.Acc_null(:,1));
            Acc_lo = nan(1,n_sub_perm);
            Acc_hi = nan(1,n_sub_perm);
            for isp = 1:n_sub_perm
                acc_perm = sort(Ca_perm.Acc_null(isp,:)); %sort(1 - Ca_perm.Loss_null);
                pval_thresh = 0.05/n_subj;
                acc_thresh_lo = acc_perm(round((pval_thresh/2) * length(acc_perm)));
                acc_thresh_hi = acc_perm(round((1 - pval_thresh/2) * length(acc_perm)));
                Acc_lo(isp) = acc_thresh_lo;
                Acc_hi(isp) = acc_thresh_hi;
                fprintf('[*] Accuracy threshold (n=%i permutations): %.4f - %.4f\n',length(acc_perm),acc_thresh_lo,acc_thresh_hi);
            end
            
        %     for i = 1:length(SubjectsSleep)
        %         sid_int = Ca.Sub_i(i);
        %     end
            x = 1:length(Ca.Acc_mean(1,:));

            %er = errorbar(x,Ca.Acc_mean(1,:),Ca.Acc_std(1,:),Ca.Acc_std(1,:));
            %bw = br(1).BarWidth;
            
            cmap = gray(length(Ca.Acc_mean(:,1))+4);
            cmap = cmap(4:(end-1),:);
            cmap = flipud(cmap);
            
            if (trig_plot_coh_only)
                Iall = 2; %length(Ca.Acc_mean(:,1));
                boffset = 0.28;
                bwidth = 0.6;
                if (state_24h)
                    pos(3) = pos(3) - 150;
                else
                    pos(3) = pos(3) - 250;
                end
            else
                Iall = 1:length(Ca.Acc_mean(:,1));
                boffset = 0.28;
                bwidth = 0.2;
            end
            
            h = figure('Position',pos,'visible','off'); hold all;
            
            % plot accuracy
            for i = Iall %1:length(Ca.Acc_mean(:,1))
                br = bar(x+(i-1)*boffset,Ca.Acc_mean(i,:),bwidth); hold on;
                br.EdgeColor = 'none';
                br.FaceColor = cmap(i,:);
                if (trig_plot_coh_only)
                    xtickvals = br.XData;
                else
                    if (i == 2)
                        xtickvals = br.XData;
                    end
                end
            end
            % plot error bars
            for i = Iall %1:length(Ca.Acc_mean(:,1))
                er = errorbar(x+(i-1)*boffset,Ca.Acc_mean(i,:),Ca.Acc_std(i,:),Ca.Acc_std(i,:)); hold on
                er.Color = [0 0 0];
                er.LineStyle = 'none';
                er.CapSize = bwidth*20;
            end
            % plot statistical thresholds
            cmap = inferno(10);
            for i = Iall
                for j = 1:length(Acc_hi)
                    hold on;
                    x_t = linspace(x(j)-0.02,x(j)+0.58,2);
                    plot(x_t,Acc_hi(j)*ones(1,length(x_t)),'-','Color',cmap(4,:));
                end
            end
            
            % plot average accuracy for both
            if (trig_plot_coh_only)
                acc_avg = mean(Ca.Acc_mean(2,:));
            else
                acc_avg = mean(Ca.Acc_mean(3,:));
            end
            ax = gca;
            plot(ax.XLim,acc_avg*[1 1],':','Color',0.5*[1 1 1]);
            %return

            set(gca, 'Layer', 'top');
            xticks(xtickvals);
            xticklabels(Ca.Sub_i);
            xlabel('Subject');
            ylabel('Accuracy');
            box off;
            set(gca,'TickDir','out');
            ylim([0.5 1]);
            if (~trig_plot_coh_only)
                legend({'Small-World','Coherence','Both'},'Location','EastOutside');
            end
            set(gcf,'renderer','Painters');
            if (state_24h)
                outfn = sprintf('./figures/T21/accuracy_conly-%i_%s_accmu-%i_metric-%i_%s',trig_plot_coh_only,Btxt{ib},round(1000*acc_avg),metrici,btxt);
            else
                outfn = sprintf('./figures/T21/accuracy_conly-%i_%s_accmu-%i_metric-%i',trig_plot_coh_only,Btxt{ib},round(1000*acc_avg),metrici);
            end
            print(h,outfn,'-dpng');
            print(h,outfn,'-depsc');
            close(h);
            
            % display in txt
            if (isnan(Beh_txt))
                Beh_txt = 'Sleep';
            end
            fprintf('%s ',Beh_txt);
            for i = Iall
                for i2 = 1:length(Ca.Acc_mean(i,:))
                    fprintf('(%i)mu:%.4f,std:%.4f ',Ca.Sub_i(i2),Ca.Acc_mean(i,i2),Ca.Acc_std(i,i2));
                end
            end
            fprintf('\n');
        end
    end
end