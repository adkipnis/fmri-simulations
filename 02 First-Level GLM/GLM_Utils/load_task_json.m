function Opts = load_task_json(Opts, Dirs)
fname = ['task-', Opts.task ,'_bold.json']; % filename in JSON extension
str = fileread([Dirs.BIDSdir, fname]); % dedicated for reading files as text
data = jsondecode(str);
Opts.PulseSequenceType =  getfield(data, 'PulseSequenceType');
Opts.TR = getfield(data, 'RepetitionTime');
Opts.ET = getfield(data, 'EchoTime');
Opts.FlipAngle = getfield(data, 'FlipAngle');
Opts.SliceTiming = getfield(data, 'SliceTiming');
end
