files = dir('*.mat');
for file = files'
    load(file.name);
    save(file.name, 'bac', 'constants', 'grid', 'init_params', 'settings', '-v7.3')
end