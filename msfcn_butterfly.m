function msfcn_fft(block)
% Level-2 MATLAB file S-Function for times two demo.
%   Copyright 1990-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $ 

  setup(block);
  
%endfunction

function setup(block)
  %% Register dialog parameter: LMS step size 
  block.NumDialogPrms = 0;
  block.DialogPrmsTunable = {};
  
  %% Register number of input and output ports
  block.NumInputPorts  = 7; % 3 complex data inputs + shift
  block.NumOutputPorts = 5; % add complex inputs

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;

  block.SampleTime = [1,0];
  
  input_port_names = {'a_re', 'a_im', 'b_re', 'b_im', 'w_re', 'w_im'};


  for k = 1:block.NumInputPorts
      block.InputPort(k).Complexity   = 'Real';
      %block.InputPort(k).DataTypeId   = 0;
      block.InputPort(k).SamplingMode = 'Sample';
      block.InputPort(k).Dimensions   = 1;
      block.InputPort(k).DirectFeedthrough = false;
  end

  for k = 1:block.NumOutputPorts
      block.OutputPort(k).Complexity   = 'Real';
      block.OutputPort(k).DataTypeId   = 0;
      block.OutputPort(k).SamplingMode = 'Sample';
      block.OutputPort(k).Dimensions   = 1;
  end

  % Set fixed-point port data types
  %%%Fix_18_17 = block.RegisterDataTypeFxpSlopeBias(1, 18, 1.0*2^-17, 0, false);
  %%%fprintf('data type id: %d\n', Fix_18_17)
  %%%block.InputPort(1).DataTypeId = Fix_18_17;
  %%%block.InputPort(2).DataTypeId = Fix_18_17;
  %%%block.InputPort(3).DataTypeId = Fix_18_17;
  %%%block.InputPort(4).DataTypeId = Fix_18_17;
  %%%block.InputPort(5).DataTypeId = Fix_18_17;
  %%%block.InputPort(6).DataTypeId = Fix_18_17;

  %%%block.OutputPort(1).DataTypeId = Fix_18_17;
  %%%block.OutputPort(2).DataTypeId = Fix_18_17;
  %%%block.OutputPort(3).DataTypeId = Fix_18_17;
  %%%block.OutputPort(4).DataTypeId = Fix_18_17;
  
  %% Set the block simStateCompliance to default (i.e., same as a built-in block)
  block.SimStateCompliance = 'DefaultSimState';

  %% Register methods
  block.RegBlockMethod('PostPropagationSetup', @DoPostPropSetup);
  block.RegBlockMethod('Outputs',                 @Outputs);
  block.RegBlockMethod('Update',  @Update);
  block.RegBlockMethod('InitializeConditions',  @InitConditions);
  %block.RegBlockMethod('SetInputPortDataType', @SetInpPortDataType);
  %block.RegBlockMethod('SetOutputPortDataType', @SetOutputPortDataType);  
  
  %% Block runs on TLC in accelerator mode.
  block.SetAccelRunOnTLC(true);
%endfunction

%function SetOutputPortDataType(block, idx, dt)
%  block.OutputPort(idx).DataTypeID = dt;

function InitConditions(block)
    N = block.NumInputPorts;
    for k=1:N
        block.Dwork(k).Data = 0;
        %block.Dwork(k).DataTypeId = 1; % TODO should be fixed-point data
        %block.Dwork(k).DataTypeId 
    end
%endfunction

function Update(block)
    N = block.NumInputPorts;
    for k=1:N
        block.Dwork(k).Data = double(block.InputPort(k).Data);
    end

function DoPostPropSetup(block)
    block.NumDworks = block.NumInputPorts;
    for k=1:block.NumDworks
        block.Dwork(k).Name = sprintf('x%d', k);
        block.Dwork(k).Dimensions = 1;
        block.Dwork(k).DatatypeID = 0;
        block.Dwork(k).Complexity = 'Real';
        block.Dwork(k).UsedAsDiscState = true;
    end
    
%endfunction

function Outputs(block)
    N = block.NumInputPorts;
    a_re = fi(block.Dwork(1).Data, 1, 18, 17);
    a_im = fi(block.Dwork(2).Data, 1, 18, 17);
    b_re = fi(block.Dwork(3).Data, 1, 18, 17);
    b_im = fi(block.Dwork(4).Data, 1, 18, 17);
    w_re = fi(block.Dwork(5).Data, 1, 18, 17);
    w_im = fi(block.Dwork(6).Data, 1, 18, 17);
    shift = 0; %TODO

    [apbw_re, apbw_im, ambw_re, ambw_im, oflow] = ...
        butterfly_fi(a_re, a_im, b_re, b_im, w_re, w_im, shift);
    
    block.OutputPort(1).Data = double(apbw_re);
    block.OutputPort(2).Data = double(apbw_im);
    block.OutputPort(3).Data = double(ambw_re);
    block.OutputPort(4).Data = double(ambw_im);
    block.OutputPort(5).Data = double(oflow);
