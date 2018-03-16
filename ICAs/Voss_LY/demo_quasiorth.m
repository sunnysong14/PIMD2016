function demo_quasiorth()
%   function QuasiOrthDemo()
%   Applies both the standard whitening algorithm and the
%   quasi-orthogonalization algorithm to noisy data.  Displays both the
%   clean and noisy versions of the data after applying the discovered
%   preprocessing matrix to all data.
%
%   In practice, the user would only have access to the noisy data.  However,
%   the noiseless data gives a better picture of what is happening.

    % Toggle to 1 in order to repeat the previous run.  Note that this saves to the
    % same file as used by the other demo script.
    UsePriorRun = 0;

    if UsePriorRun
        load('vars')
        rng(status);
    else
        status = rng;
        rng(status);
    end
    save('vars')
        

    samples = 5000
    NoiseMagnitude = 0.5
    Sigma = [1 1
            1 2] / 2.6180;
    

    S = (rand(2, samples) - .5)*2*sqrt(3);
    A = [1 1
        0 1];
    A = A;
    X = A*S;
    X2 = X + mvnrnd(zeros(2, samples)', NoiseMagnitude*Sigma)';

    figure(1);
    plot(X2(1, :), X2(2, :), 'g+',  X(1, :), X(2, :), 'k.');
    title('Raw Data (generated from uniform vs uniform)');
%     axis('square');
    axis('equal');
    legend('Noisy Data', 'Clean Data');

    [Y2, W, b] = qwhiten(X2, 'id quas-orth');
    Y = W*(X-repmat(b, 1, size(X, 2)));
    [Q, R] = qr(W*A)
    figure(2);
    plot(Y2(1, :), Y2(2, :), 'g+', Y(1, :), Y(2, :), 'k.');
    title('Data after applying quasi-orthogonalization based on the noisy data.');
%     axis('square');
    axis('equal');
    legend('Noisy Data', 'Clean Data');

    [Y2, W, b] = qwhiten(X2, 'whiten');
    [Q, R] = qr(W*A)
    Y = W*(X-repmat(b, 1, size(X, 2)));
    figure(3);
    plot(Y2(1, :), Y2(2, :), 'g+', Y(1,:), Y(2,:), 'k.');
    title('Applied standard whitening based on the noisy data.');
%     axis('square');
    axis('equal')
    legend('Noisy Data', 'Clean Data');
end