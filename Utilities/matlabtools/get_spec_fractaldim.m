function [x, y, Dfrac, N, Lmin, Lmax] = get_spec_fractaldim(data, numline, mth, polyfit_degree, varargin)

%
% GET_FRACTALDIM( N, L, CASE_TITLE, POLYFIT_DEGREE, N_PTS, PLOTFIG )
% ==================================================================
%
% Plot the N, L parameters resulting from get_NL_boxcount and deduce the fractal dimension
% Ex :
% [D_polyfit, D_num, D_slope, Lmin_polyfit, Lmax_polyfit, Lmin_num, Lmax_num] = get_fractaldim(N_cantor, L_cantor, 1, 'Cantor', 6);
%
% INPUT
% -----
%   -data           : SPEC output data, read from read_spec(filename)
%   -numline        : Line number
%   -mth            : 0: uses the mean value to get the fractal dim
%                     1: evaluates numerically the curvature and
%                        extract the linear part
%                     2: polyfits the curve and evaluate the curvature
%                        analytically to extract the linear part.
% 	-polyfit_degree	: degree of the polyfit used to approximate the fractal
%                     dim. Only used if mth=2
%
% OUTPUT
% ------
%   -x              : log10( Lmax/L )
%   -y              : log10( N )
% 	-Dfrac          : fractal dimension
%   -N              : Number of points as function of L
%	-Lmin    		: Maximal L taken into account for evaluation of Dfrac
%	-Lmax   		: Minimal L taken into account for evaluation of Dfrac
%

    % Optional input
    % ---------------
    torplane = 1;
    base_exponent = 2;
    kmin = 2;
    kmax = 8;
    newfig = 0;
    
    % Optional outputs
    % ----------------
    Lmin = 0;
    Lmax = 0;
    
    l = length(varargin);
    if mod(l,2)~=0
        error('Invalid number of input arguments')
    end
    
    for ii=1:l/2
       field = varargin{2*ii-1};
       value = varargin{2*ii  };
       
       switch field
           case 'torplane'
               torplane = value;
           case 'base_exponent'
               base_exponent = value;
           case 'kmin'
               kmin = value;
           case 'kmax'
               kmax = value;
           case 'newfig'
               newfig = value;
           otherwise
               warning(['Inknown input field: ', field, '. Ignored...'])
       end
    end
    
    
    % Read some important input
    % --------------------------
    R = data.poincare.R(numline,torplane,:);
    Z = data.poincare.Z(numline,torplane,:);
    
    s = size( R );
    n_pts = s(3);
    
    % Perform box counting
    % ---------------------
    [N, L] = get_NL_boxcount(base_exponent, kmin, kmax, R, Z);
    
    % Set up...
    idx = 1:length(N);
    
    % Get slopes between each datapoints
    xp = diff(-log(L));
    yp = diff(log(N));
    
    if newfig
       figure('Color','w')
       loglog( L, N, 'o', 'MarkerFaceColor', 'b')
       ylabel('Number of box with a point N')
       xlabel('Box size')

    end
    
    switch mth
        case 0  % mean slope

            % take the centered finite differences and having forward and backward 
            % finite diff for the last points
            xp_ = (xp(1:end-1) + xp(2:end))/2;
            yp_ = (yp(1:end-1) + yp(2:end))/2;
            slopes = [yp(1)/xp(1), yp_./xp_, yp(end)/xp(end)];
            
            % For output...
            Lmax = max(L);
            Lmin = min(L);
            
            % polyf_1d = polyfit(-log(L),log(N), 1);
            % D_slope = polyf_1d(1);
            Dfrac = mean(slopes);
            
            
        case 1  % numerical curvature
            % idea is to remove the saturation plateau
            % the -1 comes from the correspondance between the curvature and the slopes
            min_plateau_idx_num = min(idx(N>0.95*n_pts))-1;

            % take the centered finite differences and having forward and backward 
            % finite diff for the last points
            xp_ = (xp(1:end-1) + xp(2:end))/2;
            yp_ = (yp(1:end-1) + yp(2:end))/2;
            slopes = [yp(1)/xp(1), yp_./xp_, yp(end)/xp(end)];

            % curvature
            curv_num = curvature(-log(L),log(N));

            % threshold (kind of arbitrary for the moment)
            % we take the max w.r.t. a very low curvature in the case where the 
            % points are already all well aligned
            threshold_num = max(max(curv_num)/2, 0.2);

            % find min, max indices
            % there are always at least two points in between
            [idx_min_num, idx_max_num] = find_idx_min_max(curv_num, threshold_num, min_plateau_idx_num);

            % we take the best points to compute the dimension
            % the +1 are due to the correspondance between the curvature and the slopes
            % taking a polyfit yields a similar accuracy than the mean slope
            Dfrac = mean(slopes(idx_min_num+1:idx_max_num+1));

            % determine the corresponding Lmin_num and Lmax_num
            Lmin = L(idx_max_num+1);
            Lmax = L(idx_min_num+1);
            
        case 2  % polyfit curvature
    
            % 100 logarithmically spaced values exp(x) with x between [beg, end]
            Lfit = exp(linspace(log(L(1)), log(L(end)), 100));
            one_over_L = 1./Lfit;

            % the -1 comes from the correspondance between the curvature and the slopes
            min_plateau_idx_polyfit = uint64((min(idx(N>0.95*n_pts))-1)*100/length(N));

            % fit the whole data
            p = polyfit(-log(L),log(N), polyfit_degree);
            pp = polyder(p);
            ppp = polyder(pp);

            curv_polyfit = abs(polyval(ppp,log(one_over_L)))./(1+polyval(pp,log(one_over_L)).^2).^(3/2);

            % threshold (kind of arbitrary for the moment)
            threshold_polyfit = max(max(curv_polyfit)/10, 0.2);

            % find all zones (min, max indices) with curvature below threshold
            [idx_min_polyfit, idx_max_polyfit] = find_idx_min_max(curv_polyfit, threshold_polyfit, min_plateau_idx_polyfit);

            % determine the corresponding Lmin_polyfit and Lmax_polyfit associated 
            % with kmax_edge and kmin_edge respectively if idx_min_polyfit or 
            % idx_max_polyfit is empty, it is a consequence of a too low threshold !
            Lmin = 1/one_over_L(idx_max_polyfit);
            Lmax = 1/one_over_L(idx_min_polyfit);

            % taking a polyfit yields a similar accuracy than the mean slope
            Dfrac = mean(polyval(pp,log(one_over_L(idx_min_polyfit:idx_max_polyfit-1))));
            
            if newfig
               hold on;
               xx = log(Lfit);
               yy = polyval(p, -xx);
               loglog( exp(xx), exp(yy), 'LineWidth', 2 )
               ax = gca;
               yy = ax.YLim;
               loglog( Lmin_polyfit*[1,1], yy, 'k--' )
               loglog( Lmax_polyfit*[1,1], yy, 'k--' )
            end
            
        otherwise
            
            error('Invalid input')
            
    end
    
    x = log10(Lmax./L);
    y = log10( N );

end

% =========================================================================
% =========================================================================
function curv = curvature(x, y)

	%
	% CURV = CURVATURE( X, Y )
	% ========================
	%
	% Compute the curvature of the curve defined by x,y (the starting point and the endpoints have not defined curvature here, it is not checked whether the starting point and the endpoint coincide)
	%
	% INPUT
	% -----
	% 	-x,y		: coordinates of the curve
	%
	% OUTPUT
	% ------
	% 	-curv		: curvature of the curve
	%

	xp = diff(x);	xpp = diff(xp);
	yp = diff(y);	ypp = diff(yp);

	curv = zeros(1, length(x)-2);
	% take the centered finite differences
	xp_ = (xp(1:end-1) + xp(2:end))/2;
	yp_ = (yp(1:end-1) + yp(2:end))/2;

	curv = abs(xp_.*ypp - yp_.*xpp)./(xp_.^2 + yp_.^2).^(3/2);

	% can also be implemented such as the following (exact same formula actually) 
	% for i = 1:length(xp)-1
	%	% twice the area of the triangle between these 3 pts
	% 	D = norm(xp(i)*yp(i+1)-yp(i)*xp(i+1));
	% 	a = (x(i)-x(i+1))^2 + (y(i)-y(i+1))^2;
	% 	b = (x(i+1)-x(i+2))^2 + (y(i+1)-y(i+2))^2;
	% 	c = (x(i)-x(i+2))^2 + (y(i)-y(i+2))^2;
	% 	curv(i) = 2*D/(a*b*c);
	% end 
end

% ==============================================
% ==============================================
function [idx_min, idx_max] = find_idx_min_max(curv, threshold, min_plateau_idx)

	%
	% FIND_IDX_MIN_MAX( CURV, threshold, min_plateau_idx )
	% =================================================
	%
	% Find the zones (min and max indices) associated with the minimal curvature under the threshold (one of them is the best linear zone)
	%
	% INPUT
	% -----
	% 	-curv				: curvature of the curve
	% 	-threshold			: threshold for the minimal curvature
	% 	-plateau_idx		: indices of the pts near the saturation plateau at log(n_pts)
	%
	% OUTPUT
	% ------
	% 	-idx_min,idx_max	: min and max indices for the linear zone
	%

	n_curv = length(curv);

	% find all the indices such that both neighbours are below the threshold
	start_pt_2_idx_below_threshold = [];
	end_pt_2_idx_below_threshold = [];
	for i = 1:n_curv
		% finds the pts such that \_
		if i>1 && curv(i) < threshold && curv(i-1) > threshold || i==1 && curv(1) < threshold
			start_pt_2_idx_below_threshold(end+1) = i;
		end
		% finds the pts such that _/
		if i<n_curv && curv(i+1) > threshold && curv(i) < threshold || i==n_curv && curv(n_curv) < threshold
			end_pt_2_idx_below_threshold(end+1) = i;
		end
	end

	if length(end_pt_2_idx_below_threshold) > 1 
		if ~isempty(min_plateau_idx)
			% remove the plateau from the possibilities
			while start_pt_2_idx_below_threshold(end) >= min_plateau_idx && length(end_pt_2_idx_below_threshold) > 1
				start_pt_2_idx_below_threshold(end) = [];
				end_pt_2_idx_below_threshold(end) = [];
			end
		end

		% find the largest zone delta L
		[~, idx_largest_delta_L] = max(end_pt_2_idx_below_threshold-start_pt_2_idx_below_threshold);
		% as indices of curv
		idx_min = start_pt_2_idx_below_threshold(idx_largest_delta_L);
		idx_max = end_pt_2_idx_below_threshold(idx_largest_delta_L);
	else
		% as indices of curv
		idx_min = start_pt_2_idx_below_threshold;
		idx_max = end_pt_2_idx_below_threshold;
	end
end

function [N, L] = get_NL_boxcount(base_exponent, kmin, kmax, x, y)

	%
	% GET_NL_BOXCOUNT( BASE_EXPONENT, KMIN, KMAX, X, Y )
	% ==================================================
	%
	% Get the N, L parameters necessary for fractal dimension computation using box-counting
	%
	% INPUT
	% -----
	% 	-base_exponent : base of the exponential factor by which the box sizes are reduced (1 < base <= 2 is good)
	% 	-kmin	: min exponent linked to the size of the boxes (integer >= 0)
	% 	-kmax	: max exponent linked to the size of the boxes
	% 	-x,y 	: coordinates of the points in the set to analyse (y optional)
	%
	% OUTPUT
	% ------
	% 	-N		: number of box counted
	% 	-L		: size of the box counted (normalized by the max size of the boxes Lmax)
	%

	switch nargin
	case 4
		dim = 1;
	case 5
		dim = 2;
	otherwise
		error('Too many or not enough input arguments. Type "help get_NL_fractaldim" for more informations')
	end

	n_pts = length(x);
	N = zeros(1,kmax-kmin+1);
	L = zeros(1,kmax-kmin+1);
	x_min = min(x);
	x_max = max(x);
	p = 1;

	switch dim
	case 1 % fractal is represented in 1D
		Lmax = x_max-x_min;

		for k = kmin:kmax
			s = base_exponent^k;
			L(p) = Lmax/s;

			counter = zeros(1,ceil(s));

			for i = 1:n_pts
				% x-y indices of the i-th pt
				is = 1+floor((x(i)-x_min)/L(p) );

				try
					counter(1,is) = 1;
				catch
					error(';problem')
				end
			end

			N(p) = sum(counter);
			% increment p (boxsize parameter)
			p = p+1;
		end

	case 2 % fractal is represented in 2D
		y_min = min(y);
		y_max = max(y);

		% At least 1 box in each coordinate - total of 100 boxes
		% Better precision than if one takes : 
		% 	Lmax = max([x_max-x_min, y_max-y_min]);
		% However one must be careful with the memory allocation!
		Lmax_x = x_max - x_min;
		Lmax_y = y_max - y_min;
		Lmax = min([Lmax_x, Lmax_y]);
		if abs(Lmax) < 1e-7
			% due to points being possible
			Lmax = 1;
		end
		for k = kmin:kmax
			s = base_exponent^k;
			L(p) = Lmax/s;

			counter = zeros(ceil(s * Lmax_x / Lmax), ceil(s * Lmax_y / Lmax));

			for i = 1:n_pts
				% x-y indices of the i-th pt
				is = 1+floor((x(i)-x_min)/L(p) );
				js = 1+floor((y(i)-y_min)/L(p) );

				try
					counter(is,js) = 1;
				catch
					error(';problem')
				end
			end

			N(p) = sum(sum(counter));
			% increment p (boxsize parameter)
			p = p+1;
		end
    end    
end
