using julia_petra


### test Serial Comm ###

serialComm = SerialComm()

# ensure no errors or hangs
barrier(serialComm)

@test_throws AssertionError broadcast(serialComm, [1, 2, 3], 2)
@test [1, 2, 3] == broadcast(serialComm, [1, 2, 3], 1)
@test ['a', 'b', 'c'] == broadcast(serialComm, ['a', 'b', 'c'], 1)