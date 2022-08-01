%     written by: Mohamed Rabie, 2015
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>      

% PengChen: ע�ͣ��������name�ͷ���init����Ϊ������Nix_M 2020/12/02 15:52:01 
% Ŀǰ��֧�� e-n ��Ӧ���ͣ���Ҫ��ÿ��n�Ľ��������ļ�ΪLXCat��ʽ����׺Ϊtxt���洢�ڵ����ļ����ͬ�ļ��в���������txt�ļ�

% ʹ��˵��
% Xsec.energy: 1*N double�����������tot�����Ӧ������
% Xsec.tot{j}: 1*N double����j�������tot���档Xsec.tot(i)��ΪԪ����
% Xsec.exc{j}{k}: 2*K double����j�������exc��������
%                      ����Xsec.exc{j}{k}(:,1)��eV
%                      ����Xsec.exc{j}{k}(:,2)��m^2

% ע�⣺
% �������ݵ����������Ϣ������ʹ�øý���ǰֱ�Ӳ鿴�ļ�
% info���������Ӧϸ�ڣ�����Ŀ�������ڸô������Ӧ�ڿ��ط�Ӧ���ͺ������
% E_head>1e-6
classdef ImportLXCat
    
    properties
        % cell array of gas name
        name
        % cell array of directories of LXcat files
        dir
        % ionization cross sections
        ion
        % attachment cross sections
        att
        % elastic momentum transfer cross section
        ela
        % excitation cross sections
        exc
        % effective momentum transfer cross section
        eff
        % threshold energies of ionization
        ionThresh
        % threshold energies of exciations
        excThresh
        % threshold energies of exciations
        attThresh
        % total cross secion
        tot
        % vector giving discret energies
        energy
        % plot (interactive=1) or do not plot data (interactive=0)
        interactive
        % struct: more info from data file
        info
        
        target
        
    end
    
    methods

        % =====================================================================
        %> @brief initializes ImportLXCat
        %> Parameters: [dir,ion,att,ela,exc,eff,ionThresh,excThresh,tot,energy,interactive]
        % =====================================================================
        function obj = init(obj, gasName, gasDir, interactive)
            % initializes ImportLXCat by gasDir and interactive
            obj.name = gasName;
            obj.dir = gasDir;
            obj.interactive = interactive;
            E_head=1e-3; % ����������Χ
            E_end=1e4;
            % import section data
            obj = importXsections(obj);
            obj = fillThresholds(obj); % fill Eth=0 for ion/exc/att when no data
%             obj = add_key_points1(obj, E_head, E_end); %ion/exc/att����ֵ�ؼ���
            obj = getEnergy(obj); % create over all vector energy from LXcat data
            obj = prepareForFit1(obj); % addDataPoints for ion & exc & att������Emax��0��
            obj = effective2elastic(obj); %get ela from eff when ela isempty && eff isnotempty
%             obj = add_key_points2(obj, E_head, E_end); %ela����ֵ�ؼ���
            obj = prepareForFit2(obj); % addDataPoints for ela or eff������ά�ֿ�
            obj = totalXsection(obj);
            % other
            for j = 1 : length(obj.dir)
                if isempty(obj.tot)
                    warn(['LXCat Xsec ' num2str(j) ':' obj.name{j} ' are empty'])
                else
                    disp(['[INFO] import LXCat Xsec ' num2str(j) ':' obj.name{j} ' ok']);
                end
            end
            if obj.interactive == 1
                obj.plotXsections(14,'log','log');
            end
        end
        
        function obj = init_part1(obj, gasName, gasDir, interactive)
            % ֻ����ԭʼ���ݣ���������
            obj.name = gasName;
            obj.dir = gasDir;
            obj.interactive = interactive;
            % import section data
            obj = importXsections(obj);
            obj = fillThresholds(obj); % fill Eth=0 for ion/exc/att when no data
        end
        
        function obj = init_part2(obj, E_head, E_end)
            % ������ԭʼ���ݺ�����������
            obj = add_Eth_points(obj, E_head, E_end);
            obj = add_key_points1(obj, E_head, E_end); %ion/exc/att����ֵ�ؼ���
            obj = getEnergy(obj); % create over all vector energy from LXcat data
            obj = prepareForFit1(obj); % addDataPoints for ion & exc & att������Emax��0��
            % ����eff��ela����
            obj = add_key_points2(obj, E_head, E_end); %ela����ֵ�ؼ���
            obj = prepareForFit2(obj); % addDataPoints for ela or eff������ά�ֿ�
            obj = totalXsection(obj);
            if obj.interactive == 1
                obj.plotXsections(14,'log','log');
            end
        end
        % =====================================================================
        %> @brief loads all cross sections from LXcat-files
        %> Parameters: [ion,att,ela,exc,eff,ionThresh,excThresh]
        % =====================================================================
        function obj = importXsections(obj)
            % load all cross sections from LXcat-files into arrays "ion","att", "ela", "exc" and "eff" and threshold energies into arrays "ionThresh" and "excThresh"
            
            obj = ionization(obj);
            obj = attachment(obj);
            obj = elastic(obj);
            obj = excitation(obj);
            obj = effective(obj);
            
        end
        
        % =====================================================================
        %> @brief loads ionization cross sections from LXcat-files
        %> Parameters: [ion,ionThresh]
        % =====================================================================
        function obj = ionization(obj)
            % load ionization cross sections from LXcat-files into array "ion" and ionization threshold energy into "ionThresh"
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x = strfind(C{1,1}, 'IONIZATION');
                K = sum(cell2mat(x));
                
                if K == 0
                    
                    obj.ion{j} = [];
                    obj.info.ion{j}{1}='no Ionization process';
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        
                        sigma = 0;
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        
                        ind_ion = ind + 2;
                        obj.info.ion{j}{k} = get_str_process(C,ind);
                        
                        %obj.target.ion{j}{k} = get_str_target(C,ind);%%����target
                        
                        x = strfind(C{1,1}, '---');
                        i = ind;
                        while isempty(x{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(x{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.ion{j}{k} = sigma;
                        obj.ionThresh{j}{k} = str2num(cell2mat(C{1,1}(ind_ion)));
                        
                        next = last + 3 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief loads attachment cross sections from LXcat-files
        %> Parameters: [att]
        % =====================================================================
        function obj = attachment(obj)
            % load attachment cross sections from LXcat-files into array "att"
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x = strfind(C{1,1}, 'ATTACHMENT');
                xx = strfind(C{1,1}, '---');
                K = sum(cell2mat(x));
                
                if K == 0
                    
                    obj.att{j} = [];
                    obj.info.att{j}{1}='no Attachment process';
                    obj.attThresh{j}{1} =[];
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        
                        sigma = 0;
                        
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        
                        obj.info.att{j}{k} = get_str_process(C,ind);
                        
                        i = ind;
                        while isempty(xx{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(xx{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                        end
                        
                        obj.att{j}{k} = sigma;
                        obj.attThresh{j}{k} =[];
                        
                        next = last + 1 ;
                    end
                    
                end
                
            end
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief loads elastic cross sections from LXcat-files
        %> Parameters: [ela]
        % =====================================================================
        function obj = elastic(obj)
            % load elastic cross sections from LXcat-files into array "ela"
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                % ���Ȳ���������������Ϊ�ļ�ͷ�ᱻʶ��Ϊ����cell
                
                x = strfind(C{1,1}, 'ELASTIC'); %Ѱ�� ELASTIC �෴Ӧ
                K = sum(cell2mat(x));
                
                if K == 0
                    
                    obj.ela{j} = [];
                    obj.info.ela{j}{1}='no Elastic process';
                    
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        sigma = 0;
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i; % ĳ��ELASTIC �෴Ӧ��C������
                        
                        obj.info.ela{j}{k} = get_str_process(C,ind);

                        x = strfind(C{1,1}, '---'); %Ѱ��ĳ����Ӧ�����ݷ�Χָʾ
                        i = ind;
                        while isempty(x{i})
                            i = i+1;
                        end
                        first = i+1; % ĳ����Ӧ��������ʼλ����C������
                        
                        i = first; 
                        while isempty(x{i})
                            i = i+1;
                        end
                        last = i-1; % ĳ����Ӧ�����ݽ���λ����C������
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.ela{j}{k} = sigma;
                        next = last + 2 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief loads excitation cross sections from LXcat-files
        %> Parameters: [exc,excThresh]
        % =====================================================================
        function obj = excitation(obj)
            % load excitation cross sections from LXcat-files into array "exc" and excThresh threshold energy into "ionThresh"
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x  = strfind(C{1,1}, 'EXCITATION');
                xx = strfind(C{1,1}, '---');
                K = sum(cell2mat(x));
                
                if K == 0
                    
                    obj.exc{j} = [];
                    obj.info.exc{j}{1}='no Excitation process';
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        
                        sigma = 0;
                        
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        ind_exc = ind + 2;
                        
                        obj.info.exc{j}{k} = get_str_process(C,ind);
                        
                        i = ind;
                        while isempty(xx{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(xx{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.exc{j}{k} = sigma;
                        obj.excThresh{j}{k} = str2num(cell2mat(C{1,1}(ind_exc)));
                        
                        next = last + 1 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief loads effective momentum transfer cross sections from LXcat-files
        %> Parameters: [ela]
        % =====================================================================
        function obj = effective(obj)
            % loads effective cross sections from LXcat-files into array
            % "eff"
            
            for j = 1 : length(obj.dir)
                
                file = dir([obj.dir{j} '/*.txt']);
                fid = fopen([obj.dir{j} '/' file.name]);
                
                C = textscan(fid, '%s %s %s %s %s %s %s %s' );
                
                x = strfind(C{1,1}, 'EFFECTIVE');
                K = sum(cell2mat(x));
                
                if K == 0
                    x = strfind(C{1,1}, 'MOMENTUM');
                    K = sum(cell2mat(x));
                end
                
                if K == 0
                    
                    obj.eff{j} = [];
                    obj.info.eff{j}{1}='no Effective cross sections';
                else
                    
                    
                    next = 1 ;
                    for k = 1 : K
                        sigma = 0;
                        i = next ;
                        while  isempty(x{i})
                            i = i+1;
                        end
                        ind = i;
                        
                        obj.info.eff{j}{k} = get_str_process(C,ind);
                        
                        x = strfind(C{1,1}, '---');
                        i = ind;
                        while isempty(x{i})
                            i = i+1;
                        end
                        first = i+1;
                        
                        i = first;
                        while isempty(x{i})
                            i = i+1;
                        end
                        last = i-1;
                        
                        N = length(first:last);
                        for i = 1 : N
                            ind = first+i-1;
                            sigma(i,1) = str2num(cell2mat(C{1,1}(ind)));
                            sigma(i,2) = str2num(cell2mat(C{1,2}(ind)));
                            
                        end
                        
                        obj.eff{j}{k} = sigma;
                        next = last + 2 ;
                    end
                    
                end
                
            end
            
            %%
            
            fclose(fid);
            
        end
        
        % =====================================================================
        %> @brief converts effective mometum transfer cross section to
        % elastic cross section
        %> Parameters: [ela]
        % =====================================================================
        function obj = effective2elastic(obj)
            % convert effective mometum transfer to elastic cross section
            
            % gas species
            for j = 1 : length(obj.dir)
                
                if isempty(obj.ela{j}) & ~isempty(obj.eff{j})
                    
                    % collision types
                    for k = 1 : length(obj.eff{j})
                        
                        
                        energy = obj.eff{j}{k}(:,1);
                        sigma = obj.eff{j}{k}(:,2);
                        
                        % % formula from Vahedi et al (1995) p.190
                        % beta = 1/2*( energy.*log(1+energy) )./...
                        %     ( energy - log(1+energy) );
                        % % for energy = 0 --> beta = 1:
                        % beta(energy<1e-6) = 1;
                        %
                        % %sigma = beta.*sigma;
                        
                        % substract excitation cross sections
                        for l = 1 : length(obj.exc{j})
                            
                            energy_eff = obj.exc{j}{l}(:,1);
                            sigma_eff = obj.exc{j}{l}(:,2);
                            % energy_eff = [energy_eff ; 1e6];
                            % sigma_eff = [sigma_eff ; sigma_eff(end)];
                            % sigma = sigma - interp1(energy_eff,sigma_eff,energy,'linear');
                            %PengChen: ��Ϊ��� 2020/12/25 20:56:18 
                            sigma = sigma - interp1(energy_eff,sigma_eff,energy,'linear','extrap');
                            % ��ʹ��interp_linear_loglog����Ϊ�������ֵ
                        end
                        
                        % substract ionization cross sections
                        for l = 1 : length(obj.ion{j})
                            
                            energy_eff = obj.ion{j}{l}(:,1);
                            sigma_eff = obj.ion{j}{l}(:,2);
                            sigma = sigma - interp1(energy_eff,sigma_eff,energy,'linear','extrap');
                            
                        end
                        
                        obj.ela{j}{k}(:,1) = energy;
                        obj.ela{j}{k}(:,2) = max(0,sigma);
                        
                    end
                    
                end
                
            end
            
        end
        
        % =====================================================================
        %> @brief creates over all energy from LXcat data, unique and sort
        %> Parameters: [energy]
        % =====================================================================
        function obj = getEnergy(obj)
            % create over all vector energy from LXcat data
            
            % gas species
            for j = 1 : length(obj.dir)
                
                % collision types
                sigma = obj.att;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.ela;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.exc;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                sigma = obj.ion;
                for k = 1 : length(sigma{j});
                    
                    energy = sigma{j}{k}(:,1);
                    obj.energy = [obj.energy ; energy] ;
                    
                end
                
                
                % find unique values and sort them
                obj.energy = unique(obj.energy);
                obj.energy = sort(obj.energy);
                
                
            end
            
        end
        
        % =====================================================================
        %> @brief Prepare data for fit by adding data points
        %> at energy = 0 and max(energy)
        %> Parameters: [ion,att,exc]
        % =====================================================================
        function obj = prepareForFit1(obj)
            % prepare data for fit by adding data points at energy = 0 and max(energy) for arrays "ion", "att" and "exc"
            
            E_max = max(obj.energy);
            
            obj.ion  = obj.addDataPoints(obj.ion ,E_max);
            obj.att  = obj.addDataPoints(obj.att ,E_max);
            obj.exc  = obj.addDataPoints(obj.exc ,E_max);
            
        end
        
        % =====================================================================
        %> @brief prepare data for fit by adding data points
        %> at energy = 0 and max(energy)
        %> Parameters: [ela,eff]
        % =====================================================================
        function obj = prepareForFit2(obj)
            % prepare data for fit by adding data points at energy = 0 and max(energy) for arrays "ela" and "eff"
            
            E_max = max(obj.energy);
            
            if ~isempty(obj.ela)
                obj.ela  = obj.addDataPoints(obj.ela,E_max);
            end
            
            if ~isempty(obj.eff)
                obj.eff  = obj.addDataPoints(obj.eff,E_max );
            end
            
        end
        
        function obj = add_Eth_points(obj, E_head, E_end)
            % Ϊ����NIx_M��ֵ����ion/exc/attԭʼ���ݲ�����ֵ��0����
            for j = 1 : length(obj.name)
                obj.exc{j}=add_Eth_point(obj.exc{j}, obj.excThresh{j}, E_head, E_end);
                obj.ion{j}=add_Eth_point(obj.ion{j}, obj.ionThresh{j}, E_head, E_end);
                obj.att{j}=add_Eth_point(obj.att{j}, obj.attThresh{j}, E_head, E_end);
            end
        end
        
        function obj = add_key_points1(obj, E_head, E_end)
            % Ϊ����NIx_M��ֵ����ion/exc/attԭʼ���ݲ���E_head/E_end��
            for j = 1 : length(obj.name)
                obj.exc{j}=add_head_end_point(obj.exc{j}, E_head, E_end);
                obj.ion{j}=add_head_end_point(obj.ion{j}, E_head, E_end);
                obj.att{j}=add_head_end_point(obj.att{j}, E_head, E_end);
            end
        end
        
        function obj = add_key_points_ion(obj, E_head, E_end)
            % Ϊ����NIx_M��ֵ����ionԭʼ���ݲ���E_head/E_end��
            for j = 1 : length(obj.name)
                obj.ion{j}=add_head_end_point(obj.ion{j}, E_head, E_end);
            end
        end
        
        function obj = add_key_points_exc(obj, E_head, E_end)
            % Ϊ����NIx_M��ֵ����excԭʼ���ݲ���E_head/E_end��
            for j = 1 : length(obj.name)
                obj.exc{j}=add_head_end_point(obj.exc{j}, E_head, E_end);
            end
        end
        
        function obj = add_key_points2(obj, E_head, E_end)
            % Ϊ����NIx_M��ֵ����elaԭʼ���ݲ���E_head/E_end��
            for j = 1 : length(obj.name)
                obj.ela{j}=add_head_end_point(obj.ela{j}, E_head, E_end);
            end
        end
        % =====================================================================
        %> @brief adds data point to cross section sigma at energy = 0 and
        % zero-data point at energy, if sigma is empty
        % =====================================================================
        function [sigma] = addDataPoints(obj,sigma,Emax)
            % add data point to cross section sigma at energy = 0 and Emax
            
            % gas species
            for j = 1 : length(obj.dir)
                
                if isempty(sigma{j}) 
                    % �����ļ���δ�ܵ������ݣ����Ȳ�end��
                    % ��ȫ��energy��Emax������=0
                    sigma{j} = {[Emax'  0*Emax']};
                end
                
                % collision types
                for k = 1 : length(sigma{j})
                    
                    xfit = sigma{j}{k}(:,1);
                    yfit = sigma{j}{k}(:,2);
                    % ֮ǰ�Ѿ�ִ����add_head_end_point����������
                    % add addiitonal point at beginning:
                    if xfit(1) > 0 
                        % energy��һ��ֵ����0������ֵ��
                        xfit = [0 ; xfit];
                        % ����������ڵ����
                        yfit = [yfit(1) ; yfit];
                    end
                    
                    % ĩ�˲��㣬ԭ����sumXsection��
                    if Emax > xfit(end)
                        xfit = [xfit ; Emax];
                        yfit = [yfit ; yfit(end)];
                    end
                    
                    sigma{j}{k} = [xfit  yfit] ;
                    
                end
                
            end
            
        end
        
        % =====================================================================
        %> @brief  sets ionization and excitation threshold energies to zero if
        % not available
        %> Parameters: [ionThresh,excThresh]
        % =====================================================================
        function obj = fillThresholds(obj)
            % set ionization and excitation threshold energies ionThresh and excThresh to zero if not available
            
            % gas species
            for j = 1 : length(obj.dir)
                
                % ionization energies:
                if isempty(obj.ion{j})
                    obj.ionThresh{j} = {0};
                end
                
                % ionization energies:
                if isempty(obj.exc{j})
                    obj.excThresh{j} = {0};
                end
                
                if isempty(obj.att{j})
                    obj.attThresh{j} = {0};
                end
                
            end
            
        end
        
        % =====================================================================
        %> @brief  makes total cross section and corresponding energy
        % for each species
        %> Parameters: [energy,tot]
        % =====================================================================
        function obj = totalXsection(obj)
            % make total cross section "tot" and corresponding energy for each species
            
            for j = 1 : length(obj.dir)
                
                sigma_ion = sumXsection(obj,obj.ion{j},obj.energy);
                sigma_att = sumXsection(obj,obj.att{j},obj.energy);
                sigma_ela = sumXsection(obj,obj.ela{j},obj.energy);
                sigma_exc = sumXsection(obj,obj.exc{j},obj.energy);

                obj.tot{j} = sigma_ion + sigma_att + sigma_ela + sigma_exc;
                
%                 % test code
%                 figure
%                 loglog(obj.energy,sigma_ion,'-b')
%                 hold on
%                 loglog(obj.energy,sigma_att,'-.r')
%                 loglog(obj.energy,sigma_ela,'-.k')
%                 loglog(obj.energy,sigma_exc,'-.c')
%                 loglog(obj.energy,obj.tot{j},'--y')
%                 legend('ion','att','ela','exc','tot')
%                 grid on
            end
            
        end
        
        % =====================================================================
        %> @brief makes sum of one type of cross section over one species
        % =====================================================================
        function sigma_sum = sumXsection(obj,sigma,energy)
            % make sum of one type of cross section over one species
            
            if isempty(sigma)
                sigma_sum = 0;
            else
                sigma_sum = 0;
                for k = 1 : length(sigma)
                    
                    x = sigma{k}(:,1);
                    y = sigma{k}(:,2);
                    
                    sigma_add = interp1(x,y,energy);
                    if ~isnan(sigma_add)
                        sigma_sum = sigma_sum + sigma_add;
                    end
                    
                end
            end
            
        end

        % =====================================================================
        %> @brief plots fitted cross sections from LXcat
        % =====================================================================
        function plotXsections(obj,textsize,xscale,yscale)
            %plot fitted cross sections from LXcat
                % gas species
                for j = 1 : length(obj.dir)
                    
                    [fig, ax] = obj.GetFigure(xscale,yscale);
                    hold on; grid on
                    Legend = [];
                    for i = 1 : length(obj.ion{j})
                        plot(obj.ion{j}{i}(:,1),obj.ion{j}{i}(:,2),'r.-','Linewidth',2);
                        obj.text2function( obj.ion{j}{i}(:,1),obj.ion{j}{i}(:,2) , [ num2str(obj.ionThresh{j}{i}) ' eV'] , 12 )
                        Legend{end+1} = ['ion ' num2str(i)  ' Eth=' num2str(obj.ionThresh{j}{i}) ' eV'];
                    end
                    
                    for i = 1 : length(obj.att{j})
                        plot(obj.att{j}{i}(:,1),obj.att{j}{i}(:,2),'b.-','Linewidth',2);
                        obj.text2function( obj.att{j}{i}(:,1),obj.att{j}{i}(:,2) , [ num2str(obj.attThresh{j}{i}) ' eV'] , 12 )
                        Legend{end+1} = ['att ' num2str(i)  ' Eth=' num2str(obj.attThresh{j}{i}) ' eV'];
                    end
                    
                    for i = 1 : length(obj.exc{j})
                        plot(obj.exc{j}{i}(:,1),obj.exc{j}{i}(:,2),'black.-','Linewidth',2);
                        obj.text2function( obj.exc{j}{i}(:,1),obj.exc{j}{i}(:,2) , [ num2str(obj.excThresh{j}{i}) ' eV'] , 12 )
                        Legend{end+1} =  [ 'exc ' num2str(i)  ' Eth=' num2str(obj.excThresh{j}{i}) ' eV'];
                    end
                    
                    for i = 1 : length(obj.ela{j})
                        plot(obj.ela{j}{i}(:,1),obj.ela{j}{i}(:,2),'m--','Linewidth',3);
                        Legend{end+1} = 'elastic';
                    end
                    
                    for i = 1 : length(obj.eff{j})
                        if sum(obj.eff{j}{i}(:,2)) > 0
                            plot(obj.eff{j}{i}(:,1),obj.eff{j}{i}(:,2),'c.-','Linewidth',2);
                            Legend{end+1} = 'effective momentum transfer';
                        end
                    end
                    
                    plot(obj.energy,obj.tot{j},'--','Linewidth',2);
                    Legend{end+1} = 'total';
                    
                    
                    xlabel('energy [eV]','fontsize',textsize)
                    ylabel('\sigma  [m^2]','fontsize',textsize)
                    set(gca,'fontsize',textsize)
                    %title(['Xsections:' num2str(obj.dir{j})])
                    h = legend(Legend);
                    %set(h,'Location','BestOutside');
                    set(h,'Location','Best');
                    
                    axis([1e-4 1e4 1e-24 1e-18])
                    
                    title(['gas ' num2str(j) ':' obj.name{j}])
                    
                    % һ�㲻Ӧǿ�ȵȴ����continue
                    % ��Ϊ������ͬʱ�����������
%                     uicontrol('String','Continue',...
%                         'Callback','uiresume(gcbf)');
%                     uiwait(gcf);
                    disp(['Plot of cross sections gas ' num2str(j) ':' obj.name{j} ' ok']);
                    
                end
                
        end
        
        % =====================================================================
        %> @brief get figure properties
        % =====================================================================
        function [figure1, axes1] = GetFigure(obj,xscale,yscale)
            % get figure properties
            
            % set default scales: linear
            if nargin < 3
                xscale = 'linear';
                yscale = 'linear';
            end
            %widthX = 800;
            %widthY = 450;
            fontsize = 14;
            linewidth = 0.5;
            figure1 = figure('PaperSize',[29.68 20.98],...
                'WindowStyle','docked');
            %'Position', [1200-widthX 680-widthY widthX widthY]
            
            % 'PaperOrientation','landscape',...
            axes1 = axes('Parent',figure1, ...
                'XScale',xscale, 'XMinorTick','off', ...
                'YScale',yscale,'YMinorTick','on',...
                'LineWidth',linewidth, 'FontSize',fontsize);
            
            box(axes1,'on'); grid on
            xlabel('time [ns]','LineWidth',linewidth,'FontSize',fontsize);
            ylabel('amplitude [A]','LineWidth',linewidth,'FontSize',fontsize);
        end
        
        % =====================================================================
        %> @brief writes string to x- and y-data in Figure
        % =====================================================================
        function text2function(obj,x,y,string,textsize)
            % write text "string" on the maximum of function y(x)
            
            y_max = max(y);
            y_max = y_max(1);
            x_max = x(y==max(y));
            x_max = x_max(1);
            text(x_max,y_max, string ,'fontsize',textsize)
            
        end
        
        % =====================================================================
        %> @brief update_info
        % =====================================================================
        function obj = add_threshold_to_info(obj)
            % update info by adding Eth
            for j = 1 : length(obj.name)
                 for k = 1 : length(obj.exc{j})
                     obj.info.exc{j}{k}=[obj.info.exc{j}{k}...
                         ' Eth=' num2str(obj.excThresh{j}{k}) 'eV'];
                 end
                 for k = 1 : length(obj.ion{j})
                     obj.info.ion{j}{k}=[obj.info.ion{j}{k}...
                         ' Eth=' num2str(obj.ionThresh{j}{k}) 'eV'];
                 end
                 for k = 1 : length(obj.att{j})
                     obj.info.att{j}{k}=[obj.info.att{j}{k}...
                         ' Eth=' num2str(obj.attThresh{j}{k}) 'eV'];
                 end
            end
        end
        
        % =====================================================================
        %> @brief choose exc to be used
        % =====================================================================
        function obj = choose_exc(obj, j, k_used)
            % choose exc in the array k_used for the jth gas
            assert(~isempty((k_used)))
            assert(~(max(k_used)>length(obj.exc{j}))) %����k_used������������Ӧ��Ŀ
            obj.exc{j}=obj.exc{j}(k_used); 
            obj.info.exc{j}=obj.info.exc{j}(k_used);
            obj.excThresh{j}=obj.excThresh{j}(k_used); 
        end
                        
        % =====================================================================
        %> @brief output info
        % =====================================================================
        function output_info(obj)
            % output info
            for j = 1 : length(obj.name)
                disp(obj.name{j})
                disp('Elastic')
                output_info(obj.info.ela{j})
                disp('Excitation')
                output_info(obj.info.exc{j})
                disp('Ionization')
                output_info(obj.info.ion{j})
                disp('Attachment ')
                output_info(obj.info.att{j})
                disp('Effective')
                output_info(obj.info.eff{j})
                fprintf('\n')
            end
        end
        
    end
    
end

%% aid function
function output_info(in_cell)
% ���ĳ�෴Ӧ�����й��̵�info
len=length(in_cell);
if strcmp(in_cell{1}(1:2),'no')
    disp('num of process: 0')
else
    fprintf('num of process: %d\n',len)
    for k= 1: len
        disp([num2str(k) ' ' in_cell{k}])
    end
end
end

function str_process=get_str_process(C,ind)
% Ѱ��LXCat�ļ��еĹ�����Ϣ
x1 = strfind(C{1,1}, 'PROCESS'); %Ѱ��PROCESS��Ϣ
x2 = strfind(C{1,1}, 'PARAM');
i = ind;
while isempty(x1{i})
    i = i+1;
end
first_line=i; % �õ�PROCESS��Ϣ��C������
while isempty(x2{i})
    i = i+1;
end
last_line=i-1;
str_process='';
for i_line=first_line:last_line
    for i_cell_column=1:length(C)
        object_from_cell=C{1,i_cell_column}{i_line};
        if isa(object_from_cell, 'char')
            str_process=[str_process object_from_cell];
        elseif isa(object_from_cell, 'numeric')
            str_process=[str_process num2str(object_from_cell)];
        end
    end
end
end

function str_target=get_str_target(C,ind)
% Ѱ��LXCat�ļ��еĹ�����Ϣ
x1 = strfind(C{1,1}, 'TARGET'); %Ѱ��PROCESS��Ϣ
x2 = strfind(C{1,1}, 'TARGET2');
i = ind;
while isempty(x1{i})
    i = i+1;
end
first_line=i; % �õ�PROCESS��Ϣ��C������
while isempty(x2{i})
    i = i+1;
end
last_line=i-1;
str_target='';
for i_line=first_line:last_line
    for i_cell_column=1:length(C)
        object_from_cell=C{1,i_cell_column}{i_line};
        if isa(object_from_cell, 'char')
            str_target=[str_target object_from_cell];
        elseif isa(object_from_cell, 'numeric')
            str_target=[str_target num2str(object_from_cell)];
        end
    end
end
end


function [K_xsec]=add_Eth_point(K_xsec, K_thresh, E_head, E_end)
% ��ion/exc/att���й��̵���ֵ�����
if ~isempty(K_xsec)
    for k = 1 : length(K_xsec)
        % ���봦��
        energy=K_xsec{k}(:,1);
        xsec=K_xsec{k}(:,2);
        Eth=K_thresh{k};
        assert(Eth>0)
        assert(Eth>=E_head)
        assert(Eth<=E_end)
        % % ����ֵ���ݵ� 
        isequal_Eth=@(E_one) isequal_with_error(E_one,Eth,'AbsTol',1e-7);
%         if ~any(energy==Eth) % 13.6ʱ�������´���
        if ~any(arrayfun(isequal_Eth,energy)) 
            % �����ֵ�㴦�����
            energy=[energy; Eth];
            xsec=[xsec; 0];
        end
        % ����
        [energy,idx_new] = sort(energy);
        xsec=xsec(idx_new);
        % ��ֵ��֮ǰ��Ϊ0
         xsec(energy<=Eth)=0;
         % �������
         K_xsec{k}=[energy xsec];
    end
end
end

function [K_xsec]=add_head_end_point(K_xsec, E_head, E_end)
% ��������ֵ�Ŀ����������ݵ�
if ~isempty(K_xsec)
    for k = 1 : length(K_xsec)
        % ���봦��
        energy=K_xsec{k}(:,1);
        xsec=K_xsec{k}(:,2);
        assert(all(sort(unique(energy))==energy)) %��֤��������ȥ��
        assert(~any(energy<0))
        assert(~any(xsec<0))
        assert(E_head>0)
        assert(E_end>E_head)
        %%%%%%%%%%% ��
        isequal_E_head=@(E_one) isequal_with_error(E_one,E_head,'AbsTol',1e-7);
%         if ~any(energy==E_head) % �����������´���
        if ~any(arrayfun(isequal_E_head,energy)) 
             % ��head�������ݵ����=0����head�����=0
             % ��head�������ݵ����>0�����ֵ��head��
             xsec_head=interp_linear_loglog(energy,xsec,E_head);
             energy=[E_head; energy];
             xsec=[xsec_head; xsec];
         end
         % ����
         [energy,idx_new] = sort(energy);
         xsec=xsec(idx_new);
         %%%%%%%%%%% ĩ
         isequal_E_end=@(E_one) isequal_with_error(E_one,E_end,'AbsTol',1e-7);
%         if ~any(energy==E_end) % �����������´���
        if ~any(arrayfun(isequal_E_end,energy)) 
             % ��end�������ݵ����=0����end�����=0
             % ��end�������ݵ����>0�����ֵ��end��
             xsec_end=interp_linear_loglog(energy,xsec,E_end);
             energy=[energy; E_end];
             xsec=[xsec; xsec_end];
         end
         % ����
         [energy,idx_new] = sort(energy);
         xsec=xsec(idx_new);
         % �������
         K_xsec{k}=[energy xsec];
    end
end
end

% ����Ҫ����ʹ�ø���ʱ����Nix_M�е�interp_linear_loglog��������
% function [new_y]=interp_linear_loglog(known_x,known_y,new_x)