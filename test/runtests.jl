using MINE
using Test

@testset "Functions" begin
    x = collect(LinRange(0, 1, 1000));
    y = sin.(10 * pi * x) + x;
    n = size(x)[1];
    score = MINE.mine_compute_score(MINE.MINEProblem(x, y), MINE.MINEParameter(0.6, 15, 0))
    # score = MINE.mine_compute_score(x, y)

    @test MINE.mine_check_parameter(MINE.MINEParameter(0.6, 15, 0)) == 0
    @test MINE.mine_mic(score) ≈ 1.0
    @test MINE.mine_mas(score) ≈ 0.726071574374
    @test MINE.mine_mev(score) ≈ 1.0
    @test MINE.mine_mcn(score, 0) ≈ 4.58496250072
    @test MINE.mine_mcn_general(score) ≈ 4.58496250072
    @test MINE.mine_gmic(score) ≈ 0.779360251901
    @test MINE.mine_tic(score, false) ≈ 67.6612295532
end



# @testset "Convenience functions" begin
#     # build the X matrix, 8 variables, 320 samples
#     X = rand(8, 320)
#     # build the Y matrix, 4 variables, 320 samples
#     Y = rand(4, 320)

#     # compute pairwise statistics MIC_e and normalized TIC_e between samples in X,
#     # B=9, c=5
#     mic_p, tic_p =  MINE.pstats(X; alpha=9, c=5, est="mic_e")

#     # compute statistics between each pair of samples in X and Y
#     cstats =  MINE.cstats(X, Y, alpha=9, c=5, est="mic_e")
# end
