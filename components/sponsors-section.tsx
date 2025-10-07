'use client'

import { useState, useEffect } from 'react'
import { ExternalLink, Building2 } from 'lucide-react'
import Image from 'next/image'

export interface Sponsor {
  id: string
  name: string
  logo: string
  website: string
  description: string
  tier: 'platinum' | 'gold' | 'silver'
}

export default function SponsorsSection() {
  const [hoveredSponsor, setHoveredSponsor] = useState<string | null>(null)
  const [sponsors, setSponsors] = useState<Sponsor[]>([])
  const [loading, setLoading] = useState(true)

  useEffect(() => {
    const fetchSponsors = async () => {
      try {
        const response = await fetch('/api/sponsors')
        const data = await response.json()
        setSponsors(data)
      } catch (error) {
        console.error('Error fetching sponsors:', error)
      } finally {
        setLoading(false)
      }
    }

    fetchSponsors()
  }, [])

  return (
    <section className="py-24 px-4 bg-white">
      <div className="container mx-auto">
        <div className="text-center mb-16">
          <div className="flex justify-center mb-6">
            <div className="flex h-16 w-16 items-center justify-center rounded-2xl bg-gradient-to-br from-blue-600 to-purple-600">
              <Building2 className="h-8 w-8 text-white" />
            </div>
          </div>
          <h2 className="text-4xl md:text-5xl font-bold mb-4 bg-gradient-to-r from-blue-600 to-purple-600 bg-clip-text text-transparent">
            赞助伙伴
          </h2>
          <p className="text-xl text-gray-600 max-w-3xl mx-auto">
            感谢这些优秀企业对上海生物信息学中心的支持与信任
          </p>
        </div>

        {/* 照片墙效果 */}
        {loading ? (
          <div className="flex justify-center items-center py-12">
            <div className="animate-spin rounded-full h-12 w-12 border-b-2 border-blue-600"></div>
          </div>
        ) : (
          <div className="grid grid-cols-2 md:grid-cols-3 lg:grid-cols-4 xl:grid-cols-6 gap-6">
            {sponsors.map((sponsor) => (
              <div
                key={sponsor.id}
                className="group relative aspect-square bg-gray-50 rounded-2xl overflow-hidden hover:shadow-lg transition-all duration-300 hover:scale-105 cursor-pointer"
                onMouseEnter={() => setHoveredSponsor(sponsor.id)}
                onMouseLeave={() => setHoveredSponsor(null)}
                onClick={() => window.open(sponsor.website, '_blank')}
              >
                {/* 公司Logo */}
                <div className="w-full h-full flex items-center justify-center p-4">
                  <Image
                    src={sponsor.logo}
                    alt={sponsor.name}
                    width={120}
                    height={120}
                    className="w-full h-full object-contain filter grayscale group-hover:grayscale-0 transition-all duration-300"
                  />
                </div>

                {/* 悬浮时的遮罩效果 */}
                <div className="absolute inset-0 bg-gradient-to-t from-black/60 via-black/20 to-transparent opacity-0 group-hover:opacity-100 transition-all duration-300" />

                {/* 悬浮时显示的公司名称 */}
                <div className="absolute bottom-0 left-0 right-0 p-4 text-white opacity-0 group-hover:opacity-100 transition-all duration-300">
                  <h3 className="text-sm font-semibold text-center mb-1">
                    {sponsor.name}
                  </h3>
                  <div className="flex items-center justify-center text-xs text-gray-300">
                    <span>点击访问</span>
                    <ExternalLink className="ml-1 h-3 w-3" />
                  </div>
                </div>

                {/* 悬浮时显示的详细描述（更大的卡片） */}
                {hoveredSponsor === sponsor.id && (
                  <div className="absolute bottom-full left-1/2 transform -translate-x-1/2 mb-2 w-64 bg-white border border-gray-200 rounded-lg shadow-xl p-4 z-20 opacity-0 group-hover:opacity-100 transition-all duration-300">
                    <h4 className="font-semibold text-gray-900 mb-2 text-sm">
                      {sponsor.name}
                    </h4>
                    <p className="text-xs text-gray-600 leading-relaxed">
                      {sponsor.description}
                    </p>
                    <div className="mt-2 text-xs text-blue-600 flex items-center">
                      <ExternalLink className="h-3 w-3 mr-1" />
                      访问官网
                    </div>
                    {/* 小箭头指向 */}
                    <div className="absolute top-full left-1/2 transform -translate-x-1/2 w-0 h-0 border-l-4 border-r-4 border-t-4 border-transparent border-t-gray-200"></div>
                  </div>
                )}
              </div>
            ))}
          </div>
        )}

        {/* 感谢文字 */}
        <div className="mt-16 text-center">
          <p className="text-lg text-gray-600 max-w-2xl mx-auto leading-relaxed">
            他们的支持让我们能够持续为生物信息学社区提供高质量的教育资源和服务。
            如果您有兴趣成为我们的赞助伙伴，欢迎与我们联系。
          </p>
        </div>
      </div>
    </section>
  )
}