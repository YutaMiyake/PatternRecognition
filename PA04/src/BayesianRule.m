clear all; 
close all;

eigen_size = 30;

load('genderdata/16_20/trPCA_01.txt');
load('genderdata/16_20/TtrPCA_01.txt');

load('genderdata/16_20/valPCA_01.txt');
load('genderdata/16_20/TvalPCA_01.txt');

load('genderdata/16_20/tsPCA_01.txt');
load('genderdata/16_20/TtsPCA_01.txt');

train_eigens = trPCA_01(1:eigen_size,:).';
train_labels = TtrPCA_01.';
train_data = [train_labels train_eigens];

test_eigens = tsPCA_01(1:eigen_size,:).';
test_labels = TtsPCA_01.';
test_data = [test_labels test_eigens];

valid_eigens = valPCA_01(1:eigen_size,:).';
valid_labels = TvalPCA_01.';
valid_data = [valid_labels valid_eigens];

new_test_data = cat(1,valid_data,test_data);

idxs1 = train_data(:,1) == 1;
idxs2 = train_data(:,1) == 2;
male = train_data(idxs1,2:eigen_size+1);
female = train_data(idxs2,2:eigen_size+1);

%male_size = size(male,1);
%female_size = size(female,1);
train_size = size(train_data,1);
valid_size = size(valid_data,1);
test_size = size(test_data,1);
new_test_size = test_size+valid_size;

m_mean = mean(male)';
f_mean = mean(female)';
m_cov = cov(male);
f_cov = cov(female);
m_cov_inv = inv(m_cov);
f_cov_inv = inv(f_cov);

W1 = m_cov_inv*(-0.5);
W2 = f_cov_inv*(-0.5);

w1 = m_cov_inv*m_mean;
w2 = f_cov_inv*f_mean;

part11 = 0.0;
part12 = 0.0;
part21 = 0.0;
part22 = 0.0;

part11 = -0.5*m_mean'*w1;
part21 = -0.5*f_mean'*w2;

part12 = -0.5*log(det(m_cov_inv));
part22 = -0.5*log(det(f_cov_inv));
w10 = part11 + part12;
w20 = part21 + part22;


g1 = 0;
g2 = 0;
label = 0;
ctr = 0;
for i=1:new_test_size
    u = new_test_data(i,2:eigen_size+1)' - m_mean;
    g1 = u'*W1*u + w1'*u + w10;
    v = new_test_data(i,2:eigen_size+1)' - f_mean;
    g2 = v'*W2*v + w2'*v + w20;
    
    if(g1 > g2)
        label = 1;
    else
        label = 2;
    end

    if(label ~= new_test_data(i,1))
        ctr = ctr + 1;
    end
end

test_size = new_test_size
correct = test_size - ctr
accuracy = 1 - ctr/new_test_size
