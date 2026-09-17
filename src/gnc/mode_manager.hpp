#pragma once
/// @file mode_manager.hpp
/// Hysteresis-based flight-mode state machine for detumble/acquisition/nominal switching.

namespace gnc {

enum class FlightMode {
    Detumble = 1,
    RateCapture = 2,
    SlewCapture = 3,
    Nominal = 4
};

struct ModeManagerConfig {
    double enter_rate_capture_rate_rad = 0.08726646259971647;  // 5 deg/s
    double enter_slew_capture_rate_rad = 0.02617993877991494;  // 1.5 deg/s
    double enter_nominal_rate_rad      = 0.013962634015954637; // 0.8 deg/s
    double exit_nominal_rate_rad       = 0.03490658503988659;  // 2 deg/s
    double required_rate_capture_time  = 20.0;
    double required_slew_capture_time  = 4.0;
    double required_nominal_time       = 8.0;
    double min_detumble_time           = 10.0;
    double exit_rate_capture_rate_rad  = 0.10471975511965977;  // 6 deg/s
    double exit_slew_capture_rate_rad  = 0.05235987755982988;  // 3 deg/s
    bool direct_detumble_to_nominal    = false;
    bool single_acquisition_mode       = false;
    double required_nominal_exit_time  = 2.0;
};

struct ModeManager {
    explicit ModeManager(const ModeManagerConfig& cfg) : cfg(cfg) {}

    void update(double rate_norm_rad_s, bool precision_sensor_valid,
                double dt, double time_s) {
        update(rate_norm_rad_s, precision_sensor_valid, true, true, true, dt, time_s);
    }

    void update(double rate_norm_rad_s, bool precision_sensor_valid,
                bool controller_ready, bool estimator_healthy,
                double dt, double time_s) {
        update(rate_norm_rad_s, precision_sensor_valid, controller_ready, estimator_healthy, true, dt, time_s);
    }

    void update(double rate_norm_rad_s, bool precision_sensor_valid,
                bool controller_ready, bool estimator_healthy,
                bool attitude_ready,
                double dt, double time_s) {
        switched_this_step = false;
        previous_mode = mode;
        time_in_mode_s += dt;

        if (mode == FlightMode::Detumble) {
            invalid_dwell_s = 0.0;
            if (cfg.direct_detumble_to_nominal) {
                const bool nominal_ready =
                    time_s >= cfg.min_detumble_time &&
                    precision_sensor_valid &&
                    controller_ready &&
                    estimator_healthy &&
                    attitude_ready &&
                    rate_norm_rad_s < cfg.enter_nominal_rate_rad;
                valid_dwell_s = nominal_ready ? (valid_dwell_s + dt) : 0.0;
                if (valid_dwell_s >= cfg.required_nominal_time) {
                    switch_to(FlightMode::Nominal, rate_norm_rad_s, time_s);
                }
            } else {
                const bool acquisition_ready =
                    time_s >= cfg.min_detumble_time &&
                    rate_norm_rad_s < cfg.enter_rate_capture_rate_rad &&
                    estimator_healthy;

                valid_dwell_s = acquisition_ready ? (valid_dwell_s + dt) : 0.0;
                if (valid_dwell_s >= cfg.required_rate_capture_time) {
                    switch_to(FlightMode::RateCapture, rate_norm_rad_s, time_s);
                }
            }
            return;
        }

        if (mode == FlightMode::RateCapture) {
            invalid_dwell_s = 0.0;
            const bool lost_recovery =
                rate_norm_rad_s > cfg.exit_rate_capture_rate_rad || !estimator_healthy;
            if (lost_recovery) {
                valid_dwell_s = 0.0;
                switch_to(FlightMode::Detumble, rate_norm_rad_s, time_s);
                return;
            }

            const bool nominal_ready =
                precision_sensor_valid &&
                controller_ready &&
                estimator_healthy &&
                attitude_ready &&
                rate_norm_rad_s < cfg.enter_nominal_rate_rad;
            valid_dwell_s = nominal_ready ? (valid_dwell_s + dt) : 0.0;
            if (cfg.single_acquisition_mode) {
                if (valid_dwell_s >= cfg.required_nominal_time) {
                    switch_to(FlightMode::Nominal, rate_norm_rad_s, time_s);
                }
                return;
            }

            const bool slew_ready =
                precision_sensor_valid &&
                controller_ready &&
                estimator_healthy &&
                rate_norm_rad_s < cfg.enter_slew_capture_rate_rad;
            valid_dwell_s = slew_ready ? (valid_dwell_s + dt) : 0.0;
            if (valid_dwell_s >= cfg.required_slew_capture_time) {
                switch_to(FlightMode::SlewCapture, rate_norm_rad_s, time_s);
            }
            return;
        }

        if (mode == FlightMode::SlewCapture) {
            invalid_dwell_s = 0.0;
            const bool lost_recovery =
                rate_norm_rad_s > cfg.exit_rate_capture_rate_rad || !estimator_healthy;
            if (lost_recovery) {
                valid_dwell_s = 0.0;
                switch_to(FlightMode::Detumble, rate_norm_rad_s, time_s);
                return;
            }

            const bool degraded_to_rate_capture =
                !precision_sensor_valid ||
                !controller_ready ||
                rate_norm_rad_s > cfg.exit_slew_capture_rate_rad;
            if (degraded_to_rate_capture) {
                valid_dwell_s = 0.0;
                switch_to(FlightMode::RateCapture, rate_norm_rad_s, time_s);
                return;
            }

            const bool nominal_ready =
                precision_sensor_valid &&
                controller_ready &&
                attitude_ready &&
                rate_norm_rad_s < cfg.enter_nominal_rate_rad;
            valid_dwell_s = nominal_ready ? (valid_dwell_s + dt) : 0.0;
            if (valid_dwell_s >= cfg.required_nominal_time) {
                switch_to(FlightMode::Nominal, rate_norm_rad_s, time_s);
            }
            return;
        }

        const bool nominal_invalid =
            rate_norm_rad_s > cfg.exit_nominal_rate_rad ||
            !precision_sensor_valid ||
            !controller_ready ||
            !estimator_healthy ||
            !attitude_ready;
        invalid_dwell_s = nominal_invalid ? (invalid_dwell_s + dt) : 0.0;
        if (nominal_invalid && invalid_dwell_s >= cfg.required_nominal_exit_time) {
            valid_dwell_s = 0.0;
            invalid_dwell_s = 0.0;
            if (cfg.direct_detumble_to_nominal) {
                switch_to(FlightMode::Detumble, rate_norm_rad_s, time_s);
            } else if (cfg.single_acquisition_mode) {
                if (!estimator_healthy || rate_norm_rad_s >= cfg.exit_rate_capture_rate_rad) {
                    switch_to(FlightMode::Detumble, rate_norm_rad_s, time_s);
                } else {
                    switch_to(FlightMode::RateCapture, rate_norm_rad_s, time_s);
                }
            } else if (!estimator_healthy || rate_norm_rad_s >= cfg.exit_rate_capture_rate_rad) {
                switch_to(FlightMode::Detumble, rate_norm_rad_s, time_s);
            } else if (precision_sensor_valid && controller_ready && rate_norm_rad_s < cfg.exit_slew_capture_rate_rad) {
                switch_to(FlightMode::SlewCapture, rate_norm_rad_s, time_s);
            } else {
                switch_to(FlightMode::RateCapture, rate_norm_rad_s, time_s);
            }
        } else if (!nominal_invalid) {
            invalid_dwell_s = 0.0;
        }
    }

    bool is_detumble() const { return mode == FlightMode::Detumble; }
    bool is_rate_capture() const { return mode == FlightMode::RateCapture; }
    bool is_slew_capture() const { return mode == FlightMode::SlewCapture; }
    bool is_transition() const { return is_rate_capture() || is_slew_capture(); }
    bool is_acquisition() const { return is_transition(); }
    bool is_fine_pointing() const { return mode == FlightMode::Nominal; }
    bool is_nominal() const { return is_fine_pointing(); }
    bool switched() const { return switched_this_step; }
    bool switched_to_rate_capture() const {
        return switched_this_step && mode == FlightMode::RateCapture;
    }
    bool switched_to_slew_capture() const {
        return switched_this_step && mode == FlightMode::SlewCapture;
    }
    bool switched_to_acquisition() const {
        return switched_to_rate_capture() || switched_to_slew_capture();
    }
    bool switched_to_nominal() const {
        return switched_this_step && mode == FlightMode::Nominal;
    }
    bool switched_to_detumble() const {
        return switched_this_step && mode == FlightMode::Detumble;
    }

    double mode_value() const { return static_cast<double>(static_cast<int>(mode)); }
    double transition_submode_value() const {
        if (mode == FlightMode::RateCapture) return 1.0;
        if (mode == FlightMode::SlewCapture) return 2.0;
        return 0.0;
    }
    double switch_event_value() const {
        if (!switched_this_step) return 0.0;
        if (mode == FlightMode::Nominal) return 2.0;
        if (mode == FlightMode::SlewCapture) return 1.5;
        if (mode == FlightMode::RateCapture) return 1.0;
        return -1.0;
    }

    ModeManagerConfig cfg;
    FlightMode mode = FlightMode::Detumble;
    FlightMode previous_mode = FlightMode::Detumble;
    double valid_dwell_s = 0.0;
    double invalid_dwell_s = 0.0;
    double time_in_mode_s = 0.0;
    double last_switch_time_s = 0.0;
    double last_switch_rate_rad_s = 0.0;
    bool switched_this_step = false;

private:
    void switch_to(FlightMode next_mode, double rate_norm_rad_s, double time_s) {
        if (mode == next_mode) return;
        previous_mode = mode;
        mode = next_mode;
        switched_this_step = true;
        time_in_mode_s = 0.0;
        last_switch_time_s = time_s;
        last_switch_rate_rad_s = rate_norm_rad_s;
    }
};

} // namespace gnc
