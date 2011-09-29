function [epsilon] = test_case(name, dims)

epsilon = ones(dims);

% epsilon(:,(dims(2)-6)/2:(dims(2)+6)/2) = 12.25;
switch name
    case 'straight'
        w = 16;
        epsilon(:,(dims(2)-w)/2:(dims(2)+w)/2) = 12.25;

    case 'L-bend'
        w = 10;
        epsilon(1:dims(1)/2,(1.3*dims(2)-w)/2:(1.3*dims(2)+w)/2) = 12.25;
        epsilon((1.3*dims(1)-w)/2:(1.3*dims(1)+w)/2, 1:dims(2)/2) = 12.25;

    case 'fiber'
        w = [10 0.6*dims(2)];
        epsilon(1:dims(1)/2,(dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25;
        epsilon(dims(1)/2:end,(dims(2)-w(2))/2:(dims(2)+w(2))/2) = 2.25;

    case 'mim'
        w = [30 4];
        epsilon(:,(dims(2)-w(1))/2:(dims(2)+w(1))/2) = 12.25;

        epsilon(end,:) = -2.0;
        epsilon(end,(dims(2)-w(2))/2:(dims(2)+w(2))/2) = 1.0;
%         epsilon(1:2,:) = -5.0;
%         epsilon(1:2,(dims(2)-w(2))/2:(dims(2)+w(2))/2) = 1.0;
%         epsilon(1:2,1:(dims(2)+w(2))/2) = 1.0;
    case 'sp'
        w = 10;
        epsilon(1:end-2,(dims(2)-w)/2:(dims(2)+w)/2) = 12.25;

        epsilon(end-1:end,:) = -3.0;
        epsilon(end-1:end,dims(2)/2:end) = 1.0;

    case 'sp-free'
        epsilon(:,:) = -4.0;
        epsilon(:,dims(2)/2:end) = 1.0;

    otherwise
        error('Invalid option for input parameter NAME.');
end

