%% Initializing workspace

   close all force ;
   clear all ;
   clc ;
   
 %% Output log
   disp( ['================================================'] ) ;
   disp( ['= Running AR1MA1-VBEM method for GRN inference ='] ) ;
   disp( ['================================================'] ) ; 
   disp(char(10)) ; disp( ['( use [Ctrl]+[C] to abort the execution )'] );
 %% Loading the data
   disp(char(10)) ;
   disp( [' # Choosing the first input for COSLIR...'] ) ;
   [ input_file,input_path ] = uigetfile( {'*.csv','Comma Separated Value (*.csv)'},'MultiSelect','off' ) ;
   file = strsplit( input_file,'.' ) ;
   data = dataset('file',[input_path, input_file],'delimiter',',','ReadObsNames',true) ;   
   input1 = double(data) ;
%    genes = get(data,'ObsNames') ; % samples = get(data(:,2:end),'VarNames') ;
%    G = size(Y,1) ;
%    N = size(Y,2) ;
   disp(char(10)) ;
   disp( [' # Choosing the second input for COSLIR...'] ) ;
   [ input_file,input_path ] = uigetfile( {'*.csv','Comma Separated Value (*.csv)'},'MultiSelect','off' ) ;
   file = strsplit( input_file,'.' ) ;
   data = dataset('file',[input_path, input_file],'delimiter',',','ReadObsNames',true) ;   
   input2 = double(data) ;
%% Specify the parameters
prompt = 'Determine the output path and name, no input implies using default values (outfile): ';
outfile = input(prompt,'s');
if isempty(outfile)
    outfile = 'output';
end
prompt = 'Determine the bootstrap number (an integer), no input implies using default values (bootnum): ';
bootnum = input(prompt);
if isempty(bootnum)
    bootnum = 50;
end
prompt = 'Determine the lambda (a number or a list of numbers), no input implies using default values (Lambda): ';
Lambda = input(prompt);
if isempty(Lambda)
    Lambda = [1e-6];
end
prompt = 'Determine the seed, no input implies using default values (seed): ';
seed = input(prompt, 's');
if isempty(seed)
    seed = 'default';
end

opts = struct();
prompt = 'Determine the eta (a number), no input implies using default values (eta): ';
opts.eta = input(prompt);
if isempty(opts.eta)
    opts.eta = 5;
end
prompt = 'Determine the epsilon (a number), no input implies using default values (epsilon): ';
opts.epsilon = input(prompt);
if isempty(opts.epsilon)
    opts.epsilon = 1e-5;
end
prompt = 'Determine whether not silent (1/0), no input implies using default values (silent): ';
opts.silent = input(prompt);
if isempty(opts.silent)
    opts.silent = 1;
end
prompt = 'Determine the max iteration (an integer), no input implies using default values (iter_max): ';
opts.iter_max = input(prompt);
if isempty(opts.iter_max)
    opts.iter_max = 20000;
end
prompt = 'Determine the number of iteration between two output (an integer), no input implies using default values (rho_output_num): ';
opts.rho_output_num = input(prompt); 
if isempty(opts.rho_output_num)
    opts.rho_output_num = 20;
end
prompt = 'Determine the number of iteration between two update of rho (an integer), no input implies using default values (rho_update_num): ';
opts.rho_update_num = input(prompt);
if isempty(opts.rho_update_num)
    opts.rho_update_num = 20;
end
% Run oracle test%%
GRN_Construct(input1, input2, Lambda, bootnum, outfile, seed, opts)