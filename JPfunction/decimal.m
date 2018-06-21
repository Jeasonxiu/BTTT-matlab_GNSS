function decimal = decimal(value);

decimal = num2str(value, '%10.2f');
decimal = str2num(decimal);