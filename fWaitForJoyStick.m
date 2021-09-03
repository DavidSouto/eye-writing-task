function fWaitForJoyStick()

        [~, ~, ~, buttons] = WinJoystickMex(0);
        while ~buttons(1) % test 'A' is pressed            
            [~, ~, ~, buttons] = WinJoystickMex(0);
        end

end