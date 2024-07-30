def int_to_roman(input_integer):
    """Convert an integer to a Roman numeral."""
    if not isinstance(input_integer, type(1)):
        raise TypeError("expected integer, got %s" % type(input_integer))
    if not 0 < input_integer < 4000:
        raise ValueError("Argument must be between 1 and 3999")
    ints = (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
    nums = ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    result = []
    for i in range(len(ints)):
        count = int(input_integer / ints[i])
        result.append(nums[i] * count)
        input_integer -= ints[i] * count
    return ''.join(result)

def roman_to_int(roman_str):
    """Convert a Roman numeral to an integer."""
    roman_map = {
        'I': 1, 'V': 5, 'X': 10, 'L': 50,
        'C': 100, 'D': 500, 'M': 1000
    }
    integer_value = 0
    prev_value = 0

    for char in reversed(roman_str.upper()):
        value = roman_map[char]
        if value < prev_value:
            integer_value -= value
        else:
            integer_value += value
        prev_value = value

    return integer_value