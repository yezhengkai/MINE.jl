using MINE
using Test
using DelimitedFiles

# Ref: https://minepy.readthedocs.io/en/latest/python.html#first-example
@testset "Functions" begin
    x = collect(LinRange(0, 1, 1000))
    y = sin.(10 * pi * x) + x
    n = size(x)[1]
    score = MINE.mine_compute_score(MINE.MINEProblem(x, y), MINE.MINEParameter(0.6, 15, 0))

    @test_throws ErrorException("alpha must be in (0.0, 1.0] or >= 4.0") MINE.mine_check_parameter(MINE.MINEParameter(2, 15, 0))
    @test_throws ErrorException("c must be > 0.0") MINE.mine_check_parameter(MINE.MINEParameter(0.6, -1, 0))
    @test_throws ErrorException("unknown estimator") MINE.mine_check_parameter(MINE.MINEParameter(0.6, 15, 2))
    @test MINE.mine_mic(score) ≈ 1.0
    @test MINE.mine_mas(score) ≈ 0.726071574374
    @test MINE.mine_mev(score) ≈ 1.0
    @test MINE.mine_mcn(score, 0) ≈ 4.58496250072
    @test MINE.mine_mcn_general(score) ≈ 4.58496250072
    @test MINE.mine_gmic(score) ≈ 0.779360251901
    @test MINE.mine_tic(score, false) ≈ 67.6612295532
end


# Ref: https://minepy.readthedocs.io/en/latest/python.html#convenience-functions-example
@testset "Convenience functions" begin
    # Read array generated by numpy
    X = readdlm("support_files/X.txt")
    Y = readdlm("support_files/Y.txt")
    expected_mic_p = readdlm("support_files/mic_p.txt")
    expected_tic_p = readdlm("support_files/tic_p.txt")
    expected_mic_c = readdlm("support_files/mic_c.txt")
    expected_tic_c = readdlm("support_files/tic_c.txt")

    mic_p, tic_p = MINE.pstats(X; alpha=9, c=5, est="mic_e")
    mic_c, tic_c = MINE.cstats(X, Y; alpha=9, c=5, est="mic_e")

    @test mic_p ≈ expected_mic_p
    @test tic_p ≈ expected_tic_p
    @test mic_c ≈ expected_mic_c
    @test tic_c ≈ expected_tic_c
end
