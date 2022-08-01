%% Main function to generate tests
function tests = test_get_v_init
% test get_v_init
tests = functiontests(localfunctions);
end

%% Test Functions
function test_maxwellian_velocity(testCase)
% test case 'Maxwellian velocity'
v_init=get_v_init( 1e6, 'Maxwellian velocity', [3,2]);
assert(isequal([3,2],size(v_init)))

vth=1e6;
type_init='Maxwellian velocity';
dimension_vec=[1e5,1];
v_init=get_v_init( vth, type_init, dimension_vec);
num=length(v_init);
h_plot1=histogram(v_init);
bin_width=h_plot1.BinWidth;
hold on
pdf_maxwellian_velocity=@(v) exp(-v.*v/(2*vth*vth))/(sqrt(2*pi)*vth);
plot_x=-5*vth:vth/5:5*vth;
plot(plot_x,num*bin_width*pdf_maxwellian_velocity(plot_x),'-r','LineWidth',3)
legend('实际分布','预期分布')
answer = questdlg('实际与预期一致？',['人工判断' type_init],'Y','N','Y');
verifyEqual(testCase,answer,'Y')
end

%% Optional file fixtures  
function setupOnce(testCase)  % do not change function name
% set a new path, for example
end

function teardownOnce(testCase)  % do not change function name
% change back to original path, for example
end

%% Optional fresh fixtures  
function setup(testCase)  % do not change function name
% open a figure, for example
figure
end

function teardown(testCase)  % do not change function name
% close figure, for example
close
end
